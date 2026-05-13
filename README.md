# BMS 270 Computational Immunology Workshop UCSF

This codebase covers an analysis of single cell RNA-Seq and VDJ data as part of the course BMS 270. All input data and results are provided through Box prior to the workshop.

## Structure

* Both single cell RNA-Seq and VDJ data were pre-processed ahead of time and run through 10x's cellranger software. Afterwards, filtered contigs outputted from cellranger were inputted into Immcantation and the output was stored in the **data** sub-directory. 
* All data were combined and pre-processed (```0_preprocessing.Rmd```) to omit doublets and annotate cell types. These processed data were stored in the **objects** directory and serves as the input for the workshop. 
* The final output files for the last exercise are store in the results directory **4_FindCandidates**. 
* HTML files of the completed exercises and additional resources can be found in the **resources** directory

```
├── CompImmunologyWorkshopUCSF
  │   ├── 1_filter_rnaseq_and_vdj_data.Rmd
  │   ├── 2_scrnaseq_analysis.Rmd
  │   ├── 3_repertoire_analysis.Rmd
  │   ├── 3b_bonus_repertoire_exercises.Rmd   (NEW — extra exercises)
  │   ├── 4_find_candidates.Rmd
  │   └── 5_publication_vdj_networks.Rmd      (NEW — publication-quality figures)
  │
  ├── 4_FindCandidates
  │   └── all_CSF_GLIPH2
  ├── data
  │   ├── input_immcantation
  │   ├── results_immcantation
  │   └── input_antigen_discovery
  ├── figures                                  (created by script 5)
  ├── objects
  └── resources
      ├── functions_workshop.R
      └── functions_pub_networks.R             (NEW — ggraph/circlize helpers)
```

### New material (2026 refresh)

* **`3b_bonus_repertoire_exercises.Rmd`** — seven additional exercises that
  extend `3_repertoire_analysis.Rmd`: CDR3 physicochemical properties
  (hydrophobicity, charge), public-vs-private clonotypes, compartment overlap
  with Morisita-Horn and Jaccard indices, V-gene-stratified spectratypes, a
  clone-level CSF-enrichment volcano, BCR isotype-switching summary, and
  expanded-vs-unexpanded gene-module scoring (effector, exhaustion, Trm).
* **`5_publication_vdj_networks.Rmd`** — rebuilds the BCR/TCR/GLIPH networks
  using `tidygraph` + `ggraph` with a colorblind-safe Okabe-Ito palette,
  proper legends, multi-panel composition via `patchwork`, V-J pairing chord
  diagrams (`circlize`), Zipf-style clone-size rank distributions, and Hill
  diversity profiles with bootstrap 95% CIs. All figures save as vector PDFs
  at journal column widths (85 / 115 / 180 mm) via `cairo_pdf`.
* **`resources/functions_pub_networks.R`** — shared helpers (palettes,
  themes, `build_repertoire_tidygraph`, `plot_repertoire_network`,
  `plot_vj_chord`, `plot_clone_rank`, `plot_hill_diversity`, `save_pub_figure`).

## Dependencies

This workshop uses **R v4.5.0** and packages installed with **Bioconductor 3.21**. Install instructions for all required R packages and software is provided on Box. 

**CRAN R packages to install:**
```{r}
# Key packages
install.packages(c("Seurat",
                   "kableExtra",
                   "alakazam",
                   "shazam",
                   "tigger",
                   "dowser",
                   "airr",
                   "gridExtra",
                   "Matrix",
                   "ggupset"),
                   repos="https://cloud.r-project.org")

# General packages to install
install.packages(c("devtools",
                   "knitr",
                   "dplyr",
                   "tidyr",
                   "tibble",
                   "stringr",
                   "reshape2",
                   "gplots",
                   "pheatmap",
                   "ggplot2",
                   "RColorBrewer",
                   "stringdist",
                   "igraph",
                   "data.table",
                   "cowplot",
                   "ggpubr",
                   "enrichplot",
                   "ggnewscale"),
                   repos="https://cloud.r-project.org")

# Additional packages for the publication-quality figures (scripts 3b and 5)
install.packages(c("tidygraph",
                   "ggraph",
                   "graphlayouts",
                   "patchwork",
                   "circlize",
                   "ggrepel",
                   "scales",
                   "ggsci"),
                   repos="https://cloud.r-project.org")
```

**Bioconductor R Packages to install:**
```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scRepertoire")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scrapper")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("Biostrings")
BiocManager::install("glmGamPoi")
BiocManager::install("dittoSeq")
BiocManager::install("ReactomePA ")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")
```

**Install Seurat version 5.2.1 to avoid errors with cell type annotations**
```{r}
library(remotes)
remotes::install_github("satijalab/seurat", ref = "v5.2.1")
```

**Github repos to install:**

* [Azimuth](https://github.com/satijalab/azimuth "Azimuth")
* [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder "DoubletFinder")
* [impactSingleCellToolkit](https://github.com/UCSF-Wilson-Lab/impactSingleCellToolkit "impactSingleCellToolkit")