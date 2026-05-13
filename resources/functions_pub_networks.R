# ===========================================================================
# Publication-quality VDJ network & repertoire visualization helpers
# ---------------------------------------------------------------------------
# These functions accompany 5_publication_vdj_networks.Rmd and provide a
# consistent, journal-ready aesthetic for BCR/TCR repertoire visualizations.
#
# Design goals:
#   * Use tidygraph + ggraph as the modern, themeable backend (instead of
#     base-R plot.igraph), so figures inherit ggplot2 theming.
#   * Single colorblind-friendly palette family (Okabe-Ito) used across panels
#     so that compartment / patient / status colors are stable.
#   * All plotting helpers return a ggplot object so they compose cleanly with
#     patchwork for multi-panel figures.
#   * Vector PDF as primary output; high-DPI PNG as fallback.
# ===========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(RColorBrewer)
  library(data.table)
})

# ---------------------------------------------------------------------------
# 1. Color palettes
# ---------------------------------------------------------------------------
# Okabe-Ito palette (colorblind-safe)
okabe_ito <- c(
  black     = "#000000",
  orange    = "#E69F00",
  skyblue   = "#56B4E9",
  green     = "#009E73",
  yellow    = "#F0E442",
  blue      = "#0072B2",
  vermilion = "#D55E00",
  pink      = "#CC79A7",
  grey      = "#999999"
)

pub_palettes <- list(
  # Compartment palette includes both "PB" and "PBMC" because different
  # workshop functions emit different label conventions.
  compartment = c(CSF = "#0072B2", PB = "#D55E00", PBMC = "#D55E00"),
  status      = c(Healthy = "#009E73", Disease = "#CC79A7"),
  patient     = c(`28` = "#0072B2", `32` = "#E69F00"),
  isotype     = c(
    IGHM  = "#56B4E9",
    IGHD  = "#999999",
    IGHG1 = "#E69F00",
    IGHG2 = "#D55E00",
    IGHG3 = "#CC79A7",
    IGHA1 = "#009E73",
    IGHA2 = "#117733",
    IGHE  = "#882255"
  ),
  tsubtype = c(
    `CD4 T cells` = "#56B4E9",
    `CD8 T cells` = "#D55E00"
  ),
  bsubtype = c(
    `B naive`        = "#56B4E9",
    `B intermediate` = "#E69F00",
    `B memory`       = "#009E73",
    Plasmablast      = "#CC79A7"
  )
)

# ---------------------------------------------------------------------------
# 2. ggplot2 publication themes
# ---------------------------------------------------------------------------
# A clean, journal-ready ggplot2 theme. Use base_size = 8 for Nature-sized
# panels, 10 for Cell, 12 for full-page figures.
theme_pub <- function(base_size = 10, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor      = element_blank(),
      panel.grid.major      = element_line(color = "grey92", linewidth = 0.25),
      panel.border          = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.ticks            = element_line(color = "black", linewidth = 0.4),
      axis.text             = element_text(color = "black", size = rel(0.9)),
      axis.title            = element_text(face = "plain", size = rel(1.0)),
      plot.title            = element_text(face = "bold", size = rel(1.05),
                                           hjust = 0, margin = margin(b = 4)),
      plot.subtitle         = element_text(size = rel(0.9), color = "grey30",
                                           margin = margin(b = 6)),
      plot.tag              = element_text(face = "bold", size = rel(1.4)),
      legend.title          = element_text(size = rel(0.9), face = "bold"),
      legend.text           = element_text(size = rel(0.85)),
      legend.key            = element_blank(),
      legend.background     = element_blank(),
      legend.position       = "right",
      strip.background      = element_rect(fill = "grey95", color = "black",
                                           linewidth = 0.4),
      strip.text            = element_text(face = "bold", size = rel(0.9))
    )
}

# Stripped-down theme for network panels (no axes, no grid, no background).
# Title is positioned over the plot ("plot" rather than "panel") so it never
# collides with neighbouring panels when composed with patchwork, and an
# explicit plot.margin reserves space for the title + subtitle on PNG export.
theme_network_pub <- function(base_size = 10, base_family = "") {
  theme_void(base_size = base_size, base_family = base_family) +
    theme(
      plot.title.position    = "plot",
      plot.caption.position  = "plot",
      plot.title             = element_text(face = "bold", size = rel(1.05),
                                            hjust = 0.5,
                                            margin = margin(t = 2, b = 4)),
      plot.subtitle          = element_text(size = rel(0.9), color = "grey30",
                                            hjust = 0.5,
                                            margin = margin(b = 8)),
      plot.tag               = element_text(face = "bold", size = rel(1.4)),
      plot.margin            = margin(t = 14, r = 8, b = 8, l = 8),
      legend.title           = element_text(size = rel(0.9), face = "bold"),
      legend.text            = element_text(size = rel(0.85)),
      legend.position        = "right",
      legend.key             = element_blank(),
      legend.background      = element_blank()
    )
}


# Normalize compartment labels — the workshop uses "PB" in repertoire data
# frames but the createRepertoireNodeDataframe() helper parses ids and emits
# "PBMC". This helper folds them into a consistent "PB" label so palettes,
# legends, and stats all align.
normalize_compartment <- function(x) {
  x <- as.character(x)
  x[x %in% c("PBMC", "pbmc")] <- "PB"
  x
}


# ---------------------------------------------------------------------------
# 3. Build a tidygraph repertoire network from the workshop's helper DFs
# ---------------------------------------------------------------------------
# Input:
#   node_df  - data frame returned by createRepertoireNodeDataframe(); must
#              contain columns: id, clone_id (or any group column),
#              clone_size, compartment, and any optional metadata columns.
#   edge_df  - data frame returned by createRepertoireEdgeDataframe() with
#              columns 'from' and 'to' that match node_df$id.
#   extra_node_meta - optional named list of metadata vectors keyed by
#              node id, mapped onto the resulting tbl_graph.
#
# Returns a tidygraph::tbl_graph object suitable for ggraph rendering.
build_repertoire_tidygraph <- function(node_df, edge_df,
                                       extra_node_meta = NULL) {
  stopifnot(all(c("id", "clone_size") %in% colnames(node_df)))
  # Drop self-loops if any (shouldn't be any from the workshop helpers, but
  # the safety net is cheap).
  edge_df <- edge_df[edge_df$from != edge_df$to, , drop = FALSE]
  # Drop dangling edges that point to non-existent nodes
  edge_df <- edge_df[edge_df$from %in% node_df$id & edge_df$to %in% node_df$id, ,
                     drop = FALSE]

  if (!is.null(extra_node_meta)) {
    for (col_name in names(extra_node_meta)) {
      vec <- extra_node_meta[[col_name]]
      node_df[[col_name]] <- unname(vec[node_df$id])
    }
  }

  # Coerce clone_size to numeric (some legacy code stores it as character)
  node_df$clone_size <- as.numeric(node_df$clone_size)

  tbl_graph(nodes = node_df, edges = edge_df, directed = FALSE)
}


# ---------------------------------------------------------------------------
# 3b. Category-constrained, non-overlapping repertoire layout
# ---------------------------------------------------------------------------
# Repertoire networks built by the workshop helpers contain one connected
# component per clone (each clone is a clique of cells). The default
# Fruchterman-Reingold layout (i) lets nodes overlap heavily and (ii) does
# nothing to spatially separate compartments. This function fixes both:
#
#   1. Partitions the connected components (clones) into three bins by their
#      compartment composition:  CSF-pure | mixed | PB-pure.
#   2. Packs the components into a 2-D grid that has those three regions as
#      side-by-side bands, so the eye instantly reads the spatial story
#      ("CSF cluster on the left, PB cluster on the right, shared clones in
#      the middle").
#   3. Within each component lays out the cells on a circle whose radius is
#      proportional to sqrt(clone_size), so larger clones get more room and
#      no two nodes ever overlap.
#
# Inputs:
#   g                - tbl_graph or igraph object
#   compartment_col  - node attribute holding compartment labels.
#   csf_label/pb_label - the strings that mark the two extreme compartments.
#   bin_thresh       - numeric in (0.5, 1]; clones with CSF fraction >= this
#                      are "CSF-pure", <= 1 - this are "PB-pure", rest are
#                      "Mixed".
#   node_radius      - target visual radius of a single node, in layout
#                      coordinates. Used to size per-clone circles.
#   gap_factor       - extra spacing multiplier between adjacent clones.
#   aspect           - canvas aspect ratio (width / height).
#
# Returns a data frame with columns x, y (one row per node, in the order of
# the node table of `g`), and clone_bin (CSF/Mixed/PB), suitable to be used
# as a manual layout via ggraph::create_layout(g, layout = "manual", ...).
layout_repertoire_grouped <- function(g,
                                      compartment_col = "compartment",
                                      csf_label       = "CSF",
                                      pb_label        = "PB",
                                      bin_thresh      = 0.85,
                                      node_radius     = 0.025,
                                      gap_factor      = 1.18,
                                      aspect          = 1.6,
                                      seed            = 1234) {
  set.seed(seed)
  ig <- if (inherits(g, "tbl_graph")) tidygraph::as.igraph(g) else g
  nodes <- igraph::as_data_frame(ig, what = "vertices")

  # ---- 1. Components (clones) ----
  comp     <- igraph::components(ig)
  nodes$.comp_id <- comp$membership
  n_total  <- nrow(nodes)

  # ---- 2. Per-component composition summary ----
  if (!(compartment_col %in% colnames(nodes))) {
    stop(sprintf("compartment_col '%s' not found in node attributes", compartment_col))
  }
  comp_df <- nodes %>%
    dplyr::group_by(.comp_id) %>%
    dplyr::summarize(
      size      = dplyr::n(),
      csf_n     = sum(.data[[compartment_col]] == csf_label),
      pb_n      = sum(.data[[compartment_col]] == pb_label),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(
      denom    = pmax(csf_n + pb_n, 1L),
      csf_frac = csf_n / denom,
      bin      = dplyr::case_when(
        csf_frac >= bin_thresh        ~ "CSF",
        csf_frac <= (1 - bin_thresh)  ~ "PB",
        TRUE                          ~ "Mixed"
      ),
      # Internal radius for the per-clone circle of cells
      r_clone  = pmax(node_radius * 1.4,
                      sqrt(size) * node_radius * 1.05),
      # Bounding-box radius incl. node radius and inter-clone gap
      r_pack   = (r_clone + node_radius) * gap_factor
    ) %>%
    # Stable ordering inside each bin: largest clones first
    dplyr::arrange(bin, dplyr::desc(size))

  # ---- 3. Pack components inside each compartment band ----
  # The canvas is [-1, 1] horizontally; we split it 1/3 each to CSF, Mixed,
  # PB. Within each band we do a simple shelf-pack: rows of clones with the
  # row height set by the largest clone in that row.
  bin_x_ranges <- list(
    CSF   = c(-1.0, -0.34),
    Mixed = c(-0.34, 0.34),
    PB    = c(0.34, 1.0)
  )

  # Vertical canvas height derived from aspect
  y_top    <-  1 / aspect
  y_bottom <- -1 / aspect

  pack_band <- function(sub_df, x_range, y_range) {
    if (nrow(sub_df) == 0) return(sub_df)
    band_w <- diff(x_range)
    band_h <- diff(y_range)

    cx_vec <- numeric(nrow(sub_df))
    cy_vec <- numeric(nrow(sub_df))

    # Shelf packing top -> down, left -> right
    cursor_x <- x_range[1]
    cursor_y <- y_range[2]
    row_h    <- 0

    for (i in seq_len(nrow(sub_df))) {
      r <- sub_df$r_pack[i]
      # Wrap to next row if we'd overflow the band horizontally
      if (cursor_x + 2 * r > x_range[2] + 1e-9) {
        cursor_x  <- x_range[1]
        cursor_y  <- cursor_y - 2 * row_h
        row_h     <- 0
      }
      cx_vec[i] <- cursor_x + r
      cy_vec[i] <- cursor_y - r
      cursor_x  <- cursor_x + 2 * r
      row_h     <- max(row_h, r)
    }

    # If we overflowed vertically, rescale band y to fit
    used_y <- y_range[2] - (cursor_y - 2 * row_h)
    if (used_y > band_h) {
      scale <- band_h / used_y
      cy_vec <- y_range[2] - (y_range[2] - cy_vec) * scale
    }

    sub_df$cx <- cx_vec
    sub_df$cy <- cy_vec
    sub_df
  }

  packed <- dplyr::bind_rows(lapply(
    c("CSF", "Mixed", "PB"),
    function(bn) {
      pack_band(
        dplyr::filter(comp_df, bin == bn),
        bin_x_ranges[[bn]],
        c(y_bottom, y_top)
      )
    }
  ))

  # ---- 4. Place nodes on a circle around each component centre ----
  pack_lookup <- packed %>% dplyr::select(.comp_id, cx, cy, r_clone, bin)
  nodes <- nodes %>%
    dplyr::left_join(pack_lookup, by = ".comp_id")

  out_x <- numeric(n_total)
  out_y <- numeric(n_total)
  comp_ids <- unique(nodes$.comp_id)
  for (cid in comp_ids) {
    idx     <- which(nodes$.comp_id == cid)
    cx      <- nodes$cx[idx[1]]
    cy      <- nodes$cy[idx[1]]
    r_clone <- nodes$r_clone[idx[1]]
    nin     <- length(idx)
    if (nin == 1) {
      out_x[idx] <- cx
      out_y[idx] <- cy
    } else {
      # Distribute cells on a circle; start at top, go clockwise. We add a
      # small random phase so identically sized clones don't visually rhyme.
      phase   <- stats::runif(1, 0, 2 * pi)
      angles  <- phase + seq(0, 2 * pi, length.out = nin + 1)[-1]
      out_x[idx] <- cx + r_clone * cos(angles)
      out_y[idx] <- cy + r_clone * sin(angles)
    }
  }

  data.frame(
    x         = out_x,
    y         = out_y,
    clone_bin = nodes$bin,
    .comp_id  = nodes$.comp_id,
    stringsAsFactors = FALSE
  )
}


# Convenience: convert the manual layout returned above into a ggraph
# layout object (so a caller can either pass it to ggraph(layout = "manual")
# or post-process freely).
to_ggraph_layout <- function(g, manual_df) {
  ggraph::create_layout(g, layout = "manual",
                        x = manual_df$x, y = manual_df$y)
}


# ---------------------------------------------------------------------------
# 4. Plot a repertoire network with ggraph (publication style)
# ---------------------------------------------------------------------------
# Arguments:
#   g           - tbl_graph returned by build_repertoire_tidygraph()
#   color_by    - node attribute name to color by (e.g. "compartment",
#                 "patient", "STATUS", "isotype").
#   color_pal   - named character vector mapping levels to colors; defaults
#                 to pub_palettes[[color_by]] when available.
#   size_range  - numeric length-2; min/max plotting size of nodes.
#   layout      - any layout supported by ggraph (e.g. "fr", "kk", "stress",
#                 "drl"). "stress" requires graphlayouts; "fr" is the
#                 reliable fallback.
#   label_top_n - integer, label the largest N clones (set to 0 to suppress).
#   label_col   - which node attribute provides labels.
#   edge_alpha  - 0-1 transparency for edges.
#   title       - plot title.
#   subtitle    - plot subtitle.
plot_repertoire_network <- function(g,
                                    color_by    = "compartment",
                                    color_pal   = NULL,
                                    size_range  = c(1.2, 5.5),
                                    layout      = "grouped",
                                    layout_args = list(),
                                    label_top_n = 0,
                                    label_col   = "clone_id",
                                    edge_alpha  = 0.25,
                                    edge_color  = "grey40",
                                    edge_width  = 0.25,
                                    band_labels = TRUE,
                                    title       = NULL,
                                    subtitle    = NULL,
                                    seed        = 1234) {
  set.seed(seed)

  if (is.null(color_pal) && !is.null(pub_palettes[[color_by]])) {
    color_pal <- pub_palettes[[color_by]]
  }

  # ---- Resolve layout ----
  band_df <- NULL
  if (identical(layout, "grouped")) {
    # Category-constrained layout that packs clones into CSF/Mixed/PB bands
    # and lays out each clone on a circle (no node overlap).
    manual_args <- c(list(g = g, seed = seed), layout_args)
    manual_df <- do.call(layout_repertoire_grouped, manual_args)
    lyt <- create_layout(g, layout = "manual",
                         x = manual_df$x, y = manual_df$y)
    # Attach bin info for optional band annotation
    lyt$clone_bin <- manual_df$clone_bin
    band_df <- data.frame(
      bin   = c("CSF",  "Mixed", "PB"),
      x_min = c(-1.0,   -0.34,    0.34),
      x_max = c(-0.34,   0.34,    1.00),
      stringsAsFactors = FALSE
    )
  } else if (identical(layout, "stress") &&
             requireNamespace("graphlayouts", quietly = TRUE)) {
    lyt <- create_layout(g, layout = "stress")
  } else if (identical(layout, "stress")) {
    message("graphlayouts not installed; falling back to 'fr'.")
    lyt <- create_layout(g, layout = "fr")
  } else {
    lyt <- create_layout(g, layout = layout)
  }

  p <- ggraph(lyt)

  # Subtle band rectangles to make the spatial-category story obvious. Drawn
  # behind edges/nodes; muted color so they don't fight the data.
  if (!is.null(band_df) && isTRUE(band_labels)) {
    y_lo <- min(lyt$y, na.rm = TRUE) - 0.05
    y_hi <- max(lyt$y, na.rm = TRUE) + 0.10
    band_df$y_min <- y_lo
    band_df$y_max <- y_hi
    band_df$fill  <- c("#0072B2", "#999999", "#D55E00")
    p <- p +
      geom_rect(data = band_df,
                aes(xmin = x_min, xmax = x_max,
                    ymin = y_min, ymax = y_max,
                    fill = bin),
                inherit.aes = FALSE,
                color = NA, alpha = 0.06) +
      ggplot2::scale_fill_manual(
        values = c(CSF = "#0072B2", Mixed = "#999999", PB = "#D55E00"),
        guide  = "none"
      ) +
      annotate("text",
               x      = (band_df$x_min + band_df$x_max) / 2,
               y      = y_hi - 0.02,
               label  = c("CSF-restricted clones",
                          "Compartment-shared",
                          "PB-restricted clones"),
               size   = 2.6, color = "grey25", fontface = "italic")
    # ggraph uses fill on nodes; reset to a new fill scale via ggnewscale
    if (requireNamespace("ggnewscale", quietly = TRUE)) {
      p <- p + ggnewscale::new_scale_fill()
    }
  }

  p <- p +
    geom_edge_link0(edge_colour = edge_color,
                    edge_width  = edge_width,
                    edge_alpha  = edge_alpha) +
    geom_node_point(aes(size = clone_size, fill = .data[[color_by]]),
                    shape = 21, color = "grey15", stroke = 0.35,
                    alpha = 0.92) +
    scale_size_continuous(range = size_range,
                          trans  = "log10",
                          breaks = c(2, 5, 10, 25, 50, 100),
                          name   = "Clone\nsize") +
    coord_fixed(clip = "off") +
    theme_network_pub() +
    labs(title = title, subtitle = subtitle, fill = color_by)

  if (!is.null(color_pal)) {
    p <- p + scale_fill_manual(values = color_pal, na.value = "grey70",
                               name = color_by)
  } else {
    p <- p + scale_fill_brewer(palette = "Set2", na.value = "grey70",
                               name = color_by)
  }

  # Optional labels for the largest clones only
  if (label_top_n > 0 && label_col %in% colnames(as_tibble(g, "nodes"))) {
    node_tbl <- as_tibble(g, "nodes")
    node_tbl$.x <- lyt$x
    node_tbl$.y <- lyt$y
    # Drop empty/NA labels first so labelling budget goes to real clones
    node_tbl <- node_tbl[!is.na(node_tbl[[label_col]]) &
                          nzchar(as.character(node_tbl[[label_col]])), ]
    node_tbl <- node_tbl[order(-as.numeric(node_tbl$clone_size)), ]
    node_tbl <- node_tbl[!duplicated(node_tbl[[label_col]]), ]
    keep_lbl <- head(node_tbl, label_top_n)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = keep_lbl,
        aes(x = .x, y = .y, label = .data[[label_col]]),
        inherit.aes = FALSE,
        size = 2.4, color = "black",
        bg.color = "white", bg.r = 0.12,
        max.overlaps = Inf, min.segment.length = 0,
        segment.size = 0.2
      )
    } else {
      p <- p + geom_text(
        data = keep_lbl,
        aes(x = .x, y = .y, label = .data[[label_col]]),
        inherit.aes = FALSE,
        size = 2.4, color = "black"
      )
    }
  }
  p
}


# ---------------------------------------------------------------------------
# 5. V-J gene pairing chord diagram with circlize
# ---------------------------------------------------------------------------
# Renders a chord diagram (a.k.a. circos plot) of V-gene <-> J-gene pairing
# frequencies for one compartment / status / patient. Saves to a PDF and
# returns the underlying counts data frame invisibly.
#
# vdj_df   - VDJ data frame with at least columns 'v_call', 'j_call', and a
#            'locus' column (TRB / IGH typically). Pre-filter the data frame
#            to the chain you want (e.g. only TRB or only IGH heavy).
# top_v    - integer, top N V genes to show (others lumped to "other").
# top_j    - integer, top N J genes to show.
# pdf_out  - optional path to save a publication-grade PDF.
# title    - plot title.
plot_vj_chord <- function(vdj_df,
                          v_col  = "v_call",
                          j_col  = "j_call",
                          top_v  = 15,
                          top_j  = 10,
                          v_color = "#0072B2",
                          j_color = "#D55E00",
                          pdf_out = NULL,
                          pdf_width  = 6,
                          pdf_height = 6,
                          title = NULL) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("circlize is required: install.packages('circlize')")
  }

  stopifnot(v_col %in% colnames(vdj_df), j_col %in% colnames(vdj_df))

  # Trim allele info if present (everything after '*')
  vc <- tstrsplit(vdj_df[[v_col]], "\\*")[[1]]
  jc <- tstrsplit(vdj_df[[j_col]], "\\*")[[1]]

  pair_df <- data.frame(V = vc, J = jc, stringsAsFactors = FALSE) %>%
    dplyr::filter(!is.na(V), !is.na(J), V != "", J != "")

  # Keep top_v and top_j most common; lump the rest into "other"
  top_vg <- names(sort(table(pair_df$V), decreasing = TRUE))[seq_len(min(top_v, length(unique(pair_df$V))))]
  top_jg <- names(sort(table(pair_df$J), decreasing = TRUE))[seq_len(min(top_j, length(unique(pair_df$J))))]

  pair_df$V <- ifelse(pair_df$V %in% top_vg, pair_df$V, "other V")
  pair_df$J <- ifelse(pair_df$J %in% top_jg, pair_df$J, "other J")

  mat_df <- pair_df %>%
    dplyr::group_by(V, J) %>%
    dplyr::summarize(n = dplyr::n(), .groups = "drop")

  v_levels <- c(top_vg, "other V")
  j_levels <- c(top_jg, "other J")
  v_levels <- v_levels[v_levels %in% mat_df$V]
  j_levels <- j_levels[j_levels %in% mat_df$J]

  # Build matrix
  mat <- matrix(0, nrow = length(v_levels), ncol = length(j_levels),
                dimnames = list(v_levels, j_levels))
  for (i in seq_len(nrow(mat_df))) {
    mat[as.character(mat_df$V[i]), as.character(mat_df$J[i])] <- mat_df$n[i]
  }

  grid_col <- c(
    setNames(colorRampPalette(c(v_color, "#9ECAE1"))(length(v_levels)), v_levels),
    setNames(colorRampPalette(c(j_color, "#FDD0A2"))(length(j_levels)), j_levels)
  )

  draw_chord <- function() {
    circlize::circos.clear()
    circlize::circos.par(
      gap.after = c(rep(1, length(v_levels) - 1), 8,
                    rep(1, length(j_levels) - 1), 8),
      start.degree = 90,
      track.margin = c(0.01, 0.01)
    )
    circlize::chordDiagram(
      mat,
      grid.col = grid_col,
      transparency = 0.25,
      annotationTrack = "grid",
      preAllocateTracks = list(track.height = 0.06),
      directional = 0,
      diffHeight  = 0.02
    )
    circlize::circos.trackPlotRegion(
      track.index = 1, panel.fun = function(x, y) {
        xlim <- circlize::get.cell.meta.data("xlim")
        ylim <- circlize::get.cell.meta.data("ylim")
        sec  <- circlize::get.cell.meta.data("sector.index")
        circlize::circos.text(
          mean(xlim), ylim[1] + 0.3, sec,
          facing = "clockwise", niceFacing = TRUE,
          adj = c(0, 0.5), cex = 0.6
        )
      }, bg.border = NA
    )
    if (!is.null(title)) {
      title(main = title, cex.main = 1.0, font.main = 2)
    }
  }

  if (!is.null(pdf_out)) {
    pdf(pdf_out, width = pdf_width, height = pdf_height, useDingbats = FALSE)
    draw_chord()
    dev.off()
  } else {
    draw_chord()
  }

  invisible(mat_df)
}


# ---------------------------------------------------------------------------
# 6. Clone-size rank distribution (Zipf-style plot)
# ---------------------------------------------------------------------------
# Plots the rank-frequency distribution of clones per group (group_col) and
# returns a ggplot object. A log-log slope close to -1 indicates a power-law
# clone-size distribution (Zipf), classic for expanded TCR repertoires.
plot_clone_rank <- function(clone_count_df,
                            group_col    = "compartment",
                            count_col    = "cell_count",
                            color_pal    = NULL,
                            min_count    = 1,
                            title        = "Clone-size rank distribution",
                            subtitle     = "Log-log distribution; slope ~ -1 indicates power-law (Zipf)",
                            show_fit_line = TRUE) {
  if (is.null(color_pal) && !is.null(pub_palettes[[group_col]])) {
    color_pal <- pub_palettes[[group_col]]
  }
  stopifnot(group_col %in% colnames(clone_count_df),
            count_col %in% colnames(clone_count_df))

  df <- clone_count_df %>%
    dplyr::filter(.data[[count_col]] >= min_count) %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::arrange(dplyr::desc(.data[[count_col]]), .by_group = TRUE) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::ungroup()

  p <- ggplot(df, aes(x = rank, y = .data[[count_col]],
                      color = .data[[group_col]])) +
    geom_point(size = 0.7, alpha = 0.85) +
    geom_line(linewidth = 0.35, alpha = 0.5) +
    scale_x_log10(labels = scales::label_number(accuracy = 1),
                  breaks = c(1, 10, 100, 1000, 10000)) +
    scale_y_log10(labels = scales::label_number(accuracy = 1),
                  breaks = c(1, 2, 5, 10, 25, 50, 100, 250)) +
    annotation_logticks(sides = "bl", size = 0.25) +
    labs(title = title, subtitle = subtitle,
         x = "Clone rank (log10)", y = "Clone size (log10)",
         color = group_col) +
    theme_pub()

  if (!is.null(color_pal)) {
    p <- p + scale_color_manual(values = color_pal, na.value = "grey50")
  }

  if (show_fit_line) {
    p <- p + geom_smooth(method = "lm", se = FALSE,
                        linewidth = 0.4, linetype = "dashed", alpha = 0.8)
  }

  p
}


# ---------------------------------------------------------------------------
# 7. Hill diversity curves (q = 0..4) with bootstrap CIs
# ---------------------------------------------------------------------------
# Computes Hill numbers D(q) over a range of q for each group, using a simple
# bootstrap to estimate 95% CIs. Returns a ggplot object.
#
# q = 0  -> Richness (number of clones)
# q = 1  -> exp(Shannon)
# q = 2  -> 1/Simpson
# Higher q emphasizes the contribution of the most abundant clones.
hill_q <- function(p, q) {
  p <- p[p > 0]
  if (length(p) == 0) return(NA_real_)
  if (q == 1) return(exp(-sum(p * log(p))))
  (sum(p^q))^(1 / (1 - q))
}

plot_hill_diversity <- function(vdj_df,
                                group_col   = "COMPARTMENT",
                                clone_col   = "clone_id",
                                cell_col    = "unique_cell_id",
                                q_range     = seq(0, 4, by = 0.2),
                                nboot       = 100,
                                color_pal   = NULL,
                                title       = "Hill diversity profile",
                                subtitle    = "Bootstrap 95% CI; q=0 richness, q=1 exp(Shannon), q=2 1/Simpson") {

  if (is.null(color_pal) && !is.null(pub_palettes[[group_col]])) {
    color_pal <- pub_palettes[[group_col]]
  }
  stopifnot(group_col %in% colnames(vdj_df),
            clone_col %in% colnames(vdj_df))

  # Collapse to one row per cell to get per-cell clone assignments
  cell_df <- unique(vdj_df[, c(cell_col, clone_col, group_col)])

  groups <- unique(cell_df[[group_col]])
  out_list <- vector("list", length(groups))
  names(out_list) <- groups

  for (g in groups) {
    sub <- cell_df[cell_df[[group_col]] == g, , drop = FALSE]
    cts <- table(sub[[clone_col]])
    cts <- cts[cts > 0]
    if (length(cts) < 2) {
      out_list[[g]] <- NULL
      next
    }

    boot_mat <- matrix(NA_real_, nrow = nboot, ncol = length(q_range))
    for (b in seq_len(nboot)) {
      idx <- sample(seq_along(cts), size = sum(cts), replace = TRUE,
                    prob = cts / sum(cts))
      p_b <- as.numeric(table(idx))
      p_b <- p_b / sum(p_b)
      boot_mat[b, ] <- vapply(q_range, function(q) hill_q(p_b, q), numeric(1))
    }
    p_point <- as.numeric(cts) / sum(cts)
    pt_vec  <- vapply(q_range, function(q) hill_q(p_point, q), numeric(1))

    out_list[[g]] <- data.frame(
      group   = g,
      q       = q_range,
      D       = pt_vec,
      D_lo    = apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      D_hi    = apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
    )
  }
  out_df <- dplyr::bind_rows(out_list)
  if (nrow(out_df) == 0) {
    stop("No groups had >= 2 clones; nothing to plot.")
  }

  p <- ggplot(out_df, aes(x = q, y = D, color = group, fill = group)) +
    geom_ribbon(aes(ymin = D_lo, ymax = D_hi), color = NA, alpha = 0.18) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 0.9) +
    scale_y_log10() +
    annotation_logticks(sides = "l", size = 0.25) +
    labs(title = title, subtitle = subtitle,
         x = "Diversity order q",
         y = expression(""^q*D~"(Hill number)"),
         color = group_col, fill = group_col) +
    theme_pub()

  if (!is.null(color_pal)) {
    p <- p + scale_color_manual(values = color_pal) +
             scale_fill_manual(values  = color_pal)
  }
  p
}


# ---------------------------------------------------------------------------
# 8. Convenience PDF / PNG exporter
# ---------------------------------------------------------------------------
# Saves a ggplot or patchwork object to a vector PDF (and optionally a
# matched 300 DPI PNG) sized appropriately for journal column widths:
#
#   "single" = ~85 mm  (Nature/Cell single column)
#   "1.5"    = ~115 mm
#   "double" = ~180 mm (full width)
#
save_pub_figure <- function(plot_obj,
                            file_basename,
                            width_class = c("single", "1.5", "double"),
                            height_mm   = 90,
                            png_dpi     = 300,
                            also_png    = FALSE,
                            title_pad_mm = 8) {
  width_class <- match.arg(width_class)
  width_mm <- c(single = 85, `1.5` = 115, double = 180)[width_class]

  # Wrap the plot so the title/subtitle have guaranteed breathing room. This
  # prevents the "title gets clipped on PNG" problem you see when patchwork
  # collects guides next to plots.
  if (inherits(plot_obj, "patchwork") || inherits(plot_obj, "gg")) {
    plot_obj <- plot_obj &
      ggplot2::theme(plot.margin = ggplot2::margin(
        t = title_pad_mm, r = 4, b = 4, l = 4, unit = "mm"
      ))
  }

  pdf_path <- paste0(file_basename, ".pdf")
  ggplot2::ggsave(pdf_path, plot = plot_obj,
                  width = width_mm, height = height_mm, units = "mm",
                  device = cairo_pdf, bg = "white")
  message("Wrote ", pdf_path)

  if (also_png) {
    png_path <- paste0(file_basename, ".png")
    ggplot2::ggsave(png_path, plot = plot_obj,
                    width = width_mm, height = height_mm, units = "mm",
                    dpi = png_dpi, bg = "white", limitsize = FALSE)
    message("Wrote ", png_path)
  }
  invisible(c(pdf = pdf_path,
              png = if (also_png) png_path else NA_character_))
}


# ---------------------------------------------------------------------------
# 9. Compose a patchwork-friendly multi-panel figure with title padding
# ---------------------------------------------------------------------------
# Helper that wraps patchwork::wrap_plots() and applies a consistent theme
# (top margin for titles, single shared legend on the right). Use this rather
# than raw patchwork operators when exporting to PNG to avoid title clipping.
compose_panels <- function(plots,
                            ncol      = 2,
                            tag_levels = "a",
                            title      = NULL,
                            subtitle   = NULL) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork is required: install.packages('patchwork')")
  }
  fig <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(
      title    = title,
      subtitle = subtitle,
      tag_levels = tag_levels,
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", size = 12,
                                              hjust = 0,
                                              margin = ggplot2::margin(b = 4)),
        plot.subtitle = ggplot2::element_text(size = 10, color = "grey30",
                                              hjust = 0,
                                              margin = ggplot2::margin(b = 8)),
        plot.tag      = ggplot2::element_text(face = "bold", size = 12)
      )
    ) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right",
                   plot.margin = ggplot2::margin(t = 12, r = 6, b = 6, l = 6))
  fig
}
