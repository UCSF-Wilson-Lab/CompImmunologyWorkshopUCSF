# BMS270 Workshop functions
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringdist))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggrepel))

# 1. Pre-Processing --------------------------------------

# createCloneCountTable - input is immcantation results
createCloneCountTable <- function(results,cellbarcode_col = "unique_cell_id") {
  clonotype_list          <- unique(results$clone_id)
  
  # Initialize Clonotype count table
  clone_count_df            <- as.data.frame(matrix(0,nrow = length(clonotype_list),ncol = 7))
  names(clone_count_df)     <- c("cell_count","CSF_count","PB_count","pct_total","pct_CSF","pct_PB","CSF_PB_ratio")
  row.names(clone_count_df) <- clonotype_list
  
  
  # Cell count totals
  results_csf <- results[results$COMPARTMENT %in% "CSF",]
  results_pb  <- results[results$COMPARTMENT %in% "PB",]
  
  total_cells       <- length(unique(results[,cellbarcode_col]))
  total_cells_csf   <- length(unique(results_csf[,cellbarcode_col])) 
  total_cells_pb    <- length(unique(results_pb[,cellbarcode_col])) 
  
  for(clone in clonotype_list){
    clone_df         <- results[results$clone_id %in% clone,]
    clone_cell_count <- length(unique(clone_df[,cellbarcode_col]))
    
    # Determine CSF and PB count
    clone_csf_count <- 0
    clone_pb_count  <- 0
    clone_df_csf <- clone_df[clone_df$COMPARTMENT %in% "CSF",]
    clone_df_pb  <- clone_df[clone_df$COMPARTMENT %in% "PB",]
    
    if(nrow(clone_df_csf) > 0){
      clone_csf_count <- length(unique(clone_df_csf[,cellbarcode_col]))
    }
    if (nrow(clone_df_pb) > 0){
      clone_pb_count <- length(unique(clone_df_pb[,cellbarcode_col]))
    }
    
    # Calculate percentage of total clones in the patient
    clone_pct     <- (clone_cell_count / total_cells) * 100
    clone_pct_csf <- (clone_csf_count / total_cells) * 100
    clone_pct_pb  <- (clone_pb_count / total_cells) * 100
    
    # Calculate CSF:PB ratio
    clone_pct_pb_ratio <- clone_pct_pb
    if(clone_pb_count == 0){
      clone_pct_pb_ratio  <- (1 / total_cells_pb) * 100
    }
    csf_pb_ratio <- clone_pct_csf / clone_pct_pb_ratio
    
    # Add clonotype results to data frame
    results_row <- clone_count_df[row.names(clone_count_df) %in% clone,]
    results_row$pct_total    <- clone_pct
    results_row$pct_CSF      <- clone_pct_csf
    results_row$pct_PB       <- clone_pct_pb
    results_row$cell_count   <- clone_cell_count
    results_row$CSF_count    <- clone_csf_count
    results_row$PB_count     <- clone_pb_count
    results_row$CSF_PB_ratio <- csf_pb_ratio
    
    clone_count_df[row.names(clone_count_df) %in% clone,] <- results_row
    
    
  }
  
  return(clone_count_df)
}

# createPatientCloneCountTable 
#  - input is immcantation results
#  - Patient level version of createCloneCountTable
createPatientCloneCountTable <- function(results,cellbarcode_col = "unique_cell_id",patient_col = "patient",timepoint_col = NULL) {
  clonotype_vec  <- results$clone_id
  patient_vec    <- results[,patient_col]
  clonotype_list <- paste(patient_vec,clonotype_vec,sep = "-")
  if(!is.null(timepoint_col)){
    clonotype_list <- paste(clonotype_list,results[,timepoint_col],sep = "-")
  }
  
  uni_patients <- tstrsplit(clonotype_list,"-")[[1]] %>% unique()
  
  clonotype_list          <- unique(clonotype_list)
  
  
  # Initialize Clonotype count table
  clone_count_df            <- as.data.frame(matrix(0,nrow = length(clonotype_list),ncol = 7))
  names(clone_count_df)     <- c("cell_count","CSF_count","PB_count","pct_total","pct_CSF","pct_PB","CSF_PB_ratio")
  row.names(clone_count_df) <- clonotype_list
  
  for (pt in uni_patients) {
    results_pt <- results[ results[,patient_col] %in% pt,]
    
    # Cell count totals
    results_csf <- results_pt[results_pt$COMPARTMENT %in% "CSF",]
    results_pb  <- results_pt[results_pt$COMPARTMENT %in% "PB",]
    
    total_cells       <- length(unique(results_pt[,cellbarcode_col]))
    total_cells_csf   <- length(unique(results_csf[,cellbarcode_col])) 
    total_cells_pb    <- length(unique(results_pb[,cellbarcode_col])) 
    
    # filter clonotypes that contain this patient
    pt_clonotype_list <- clonotype_list
    names(pt_clonotype_list) <- tstrsplit(pt_clonotype_list,"-")[[1]]
    pt_clonotype_list        <- pt_clonotype_list[names(pt_clonotype_list) %in% pt] %>% as.character()
    
    
    for(clone_str in pt_clonotype_list){
      clone            <- tstrsplit(clone_str,"-")[[2]]
      ### timepoint
      tp <- NULL
      if(!is.null(timepoint_col)){tp <- tstrsplit(clone_str,"-")[[3]]}
      
      clone_df         <- results_pt[results_pt$clone_id %in% clone,]
      if(!is.null(tp)){
        clone_df <- clone_df[ clone_df[,timepoint_col] %in% tp,]
      }
      clone_cell_count <- length(unique(clone_df[,cellbarcode_col]))
      
      # Determine CSF and PB count
      clone_csf_count <- 0
      clone_pb_count  <- 0
      clone_df_csf <- clone_df[clone_df$COMPARTMENT %in% "CSF",]
      clone_df_pb  <- clone_df[clone_df$COMPARTMENT %in% "PB",]
      
      if(nrow(clone_df_csf) > 0){
        clone_csf_count <- length(unique(clone_df_csf[,cellbarcode_col]))
      }
      if (nrow(clone_df_pb) > 0){
        clone_pb_count <- length(unique(clone_df_pb[,cellbarcode_col]))
      }
      
      # Calulate percentage of total clones in the patient
      clone_pct     <- (clone_cell_count / total_cells) * 100  # patient/timepoint level denominators
      clone_pct_csf <- (clone_csf_count / total_cells_csf) * 100
      clone_pct_pb  <- (clone_pb_count / total_cells_pb) * 100
      
      # Calculate CSF:PB ratio
      clone_pct_pb_ratio <- clone_pct_pb
      if(clone_pb_count == 0){
        clone_pct_pb_ratio  <- (1 / total_cells_pb) * 100
      }
      csf_pb_ratio <- clone_pct_csf / clone_pct_pb_ratio
      
      # Add clonotype results to data frame
      results_row <- clone_count_df[row.names(clone_count_df) %in% clone_str,]
      results_row$pct_total    <- clone_pct
      results_row$pct_CSF      <- clone_pct_csf
      results_row$pct_PB       <- clone_pct_pb
      results_row$cell_count   <- clone_cell_count
      results_row$CSF_count    <- clone_csf_count
      results_row$PB_count     <- clone_pb_count
      results_row$CSF_PB_ratio <- csf_pb_ratio
      
      clone_count_df[row.names(clone_count_df) %in% clone_str,] <- results_row
    }
  }
  
  # Add it patient and clone columns
  clone_count_df$patient <- tstrsplit(row.names(clone_count_df),"-")[[1]]
  if(!is.null(timepoint_col)){
    clone_count_df$timepoint <- tstrsplit(row.names(clone_count_df),"-")[[3]]
  }
  clone_count_df$clone_id <- tstrsplit(row.names(clone_count_df),"-")[[2]]
  
  return(clone_count_df)
}


# addExpansionCategories - For each clonotype add an expansion category
#    Highly expanded - cell count > 5
#    Moderate        - 1 < cells count < 5
#    Non-expanded    - cell count = 1
#    results column  = EXPANSION_CATEGORY 
addExpansionCategories <- function(results,clone_count_df) {
  # Transfer clonotype Freq back to full results
  clonotype_count   <- results$clone_id
  
  for(clone in row.names(clone_count_df)){
    clone_row  <- clone_count_df[clone,]
    cell_count <- clone_row$cell_count
    clonotype_count[clonotype_count %in% clone] <- cell_count
  }
  
  # Create different categories for 
  clonotype_count[clonotype_count %in% "1"]                            <- "Non-expanded"
  clonotype_count[clonotype_count %in% as.character(2:5)]              <- "Moderate"
  clonotype_count[! clonotype_count %in% c("Non-expanded","Moderate")] <- "High"
  
  results$EXPANSION_CATEGORY <- clonotype_count
  
  return(results)
}


# 3. Analysis --------------------

plotDGEvolcano <- function(results_dge, YLIM = c(-10,300),title = "DGE",p.adj.thresh = 0.05,cluster_filt = NULL,
                          with_gene_labels = TRUE,show_top_genes = 20) {
  # filter results based on comparison
  if(!is.null(cluster_filt)){results_dge <- results_dge[results_dge$cluster %in% cluster_filt,]}
  
  # Annotate Significant vs. not significant
  sig_genes <- results_dge$p_val_adj
  sig_genes[sig_genes < p.adj.thresh] <- "Significant"
  sig_genes[! sig_genes %in% "Significant"] <- "Not Significant"
  results_dge$sig_genes <- sig_genes
  
  # Scale P-adj values
  results_dge$p_val_adj_log10 <- -log10(results_dge$p_val_adj)
  results_dge$p_val_adj_log10[results_dge$p_val_adj_log10 == Inf] <- -10
  
  # Plot
  p <- ggplot(results_dge) +
    geom_point(aes(x=avg_log2FC, y=p_val_adj_log10, colour=sig_genes)) +
    scale_color_manual(values=c("black", "red")) + 
    ggtitle(title) +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) +
    scale_y_continuous(limits = YLIM) +
    scale_x_continuous(limits = c(-1.2,1.2))
  
  # Label top genes if specified
  if(with_gene_labels){
    # Subset top genes by p adj value
    results_dge_sort <- results_dge[order(results_dge$p_val_adj_log10,decreasing = T),]
    #if(!is.null(cluster_filt)){results_dge_sort <- results_dge_sort[results_dge_sort$cluster %in% cluster_filt,]}
    results_dge_top  <- results_dge_sort[1:show_top_genes,]
    p <- p+geom_text_repel(data=results_dge_top, aes(label= gene,x=avg_log2FC, y=p_val_adj_log10))
  }
  
  return(p)
}

# 4. Network Analysis Functions --------------------

### edge_list_from_distance.R
# find the clusters for each germline usage and build the edge list out of it,
# then combine all the germline together
cluster_by_distance_tcr <-function(data, distance_method, distance_cutoff){
  #CDR3<-as.character(data$CDR3)
  CDR3<-paste(data$JUNCTION_10X_AA, data$CDR3_AA_alpha,sep = ":")
  data_matrix<-stringdistmatrix(CDR3,CDR3,method= distance_method)
  
  # this is not necessary.  --Hao
  # data_matrix<-as.matrix(data_matrix) # when as.matrix is necessay? how to decide if a var is matrix?
  # TODO: distance cutoff, distance method
  
  # two cdr3 is connected if the distance less or equal to 1
  data_matrix<-ifelse(data_matrix<=distance_cutoff,1,0)
  
  ## remove self connected edge
  for (i in 1:nrow(data_matrix)){
    data_matrix[i,i]<-0
  }
  
  rownames(data_matrix)<-data$Node_ID
  colnames(data_matrix)<-data$Node_ID
  
  
  ## format the distance matrix into the adjacency matrix 
  data_adj<-graph.adjacency(data_matrix,mode=c("undirected"))
  
  ## get the edge list from the adjacency matrix
  
  data.edgelist<-as.data.frame(get.edgelist(data_adj))
  
  return(data.edgelist)
}



# Add a clone count column to VDJ results
addCloneCountColToVDJresults <- function(results, clone_counts, target.col = "", new.col.name = "") {
  conversion_vec        <- clone_counts[,target.col]
  names(conversion_vec) <- row.names(clone_counts)
  
  new_col    <- results$clone_id
  
  clone_list <- unique(new_col)
  for (clone in clone_list) {
    val <- unique(conversion_vec[clone])
    new_col[new_col %in% clone] <- val
  }
  
  results[,new.col.name] <- new_col
  return(results)
}



# Collapse mAb table to 1 line per cell (alpha and Beta)
formatTCRmabTable <- function(mab_table) {
  # Subset only keep cells with Alpha and Beta
  mab_table$v_call <- tstrsplit(mab_table$v_call,"\\*")[[1]]
  mab_table$j_call <- tstrsplit(mab_table$j_call,"\\*")[[1]]
  
  mab_table_beta  <- mab_table[mab_table$locus %in% c("TRB"),]
  row.names(mab_table_beta) <- mab_table_beta$CELL_BARCODES
  
  mab_table_alpha <- mab_table[mab_table$locus %in% c("TRA"),]
  alpha_barcodes  <- mab_table_alpha$CELL_BARCODES
  mab_table_alpha <- mab_table_alpha[,c("v_call","j_call","junction_aa")]
  names(mab_table_alpha) <- c("vgene_imm_alpha","jgene_imm_alpha","CDR3_AA_alpha")
  
  row.names(mab_table_alpha) <- alpha_barcodes
  
  
  mab_table       <- merge(mab_table_beta,mab_table_alpha, by=0,all=TRUE)
  mab_table <- mab_table[order(mab_table[,2]),]
  row.names(mab_table) <- mab_table[,1]
  mab_table[,1] <- NULL
  
  # Add Cluter Criteria Column (Alpha:Beta)
  cluster_criteria_alpha     <- paste(mab_table$vgene_imm_alpha,mab_table$jgene_imm_alpha,mab_table$CDR3_AA_alpha,sep="_")
  cluster_criteria_beta      <- paste(mab_table$v_call,mab_table$j_call,mab_table$junction_aa,sep="_")
  mab_table$cluster_criteria <- paste(cluster_criteria_alpha,cluster_criteria_beta, sep=":")
  
  # Add Node ID for network analysis
  #mab_table$Node_ID <-paste(mab_table$CLONE,mab_table$CELL_BARCODES,sep=":")
  mab_table$Node_ID  <- paste(mab_table$clone_id,mab_table[,1],sep = ":")
  # Remove NAs
  mab_table <- mab_table[!is.na(mab_table$locus),]
  
  return(mab_table)
}


# 4. Candidate Antigen Plot Functions --------------------

# peptide clustering
getGroupPeptideOrderHeatmap <- function(plot_df,select_samples,select_peptides = NULL) {
  # Get matrix of sample subset
  plot_df     <- plot_df[,names(plot_df) %in% select_samples]
  if(!is.null(select_peptides)){
    plot_df     <- plot_df[row.names(plot_df) %in% select_peptides,]
  }
  heatmap_mtx <- as.matrix(plot_df)
  
  # 1 - Heirarchical clustering - all peptides first 
  d        <- dist(heatmap_mtx)
  hc       <- hclust(d)
  hc_peptide_order <- rownames(heatmap_mtx)[hc$order]
  
  return(hc_peptide_order)
}


# sample clustering
getGroupSampleOrderHeatmap <- function(plot_df,select_samples,select_peptides = NULL) {
  # Get matrix of sample subset
  plot_df     <- plot_df[,names(plot_df) %in% select_samples]
  if(!is.null(select_peptides)){
    plot_df <- plot_df[row.names(plot_df) %in% select_peptides,]
  }
  heatmap_mtx <- t(as.matrix(plot_df))
  
  # Heirarchical clustering 
  d        <- dist(heatmap_mtx)
  hc       <- hclust(d)
  hc_sample_order <- rownames(heatmap_mtx)[hc$order]
  
  return(hc_sample_order)
}

