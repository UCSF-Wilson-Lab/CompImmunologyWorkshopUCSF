# BMS270 Workshop functions
library(rjson)     
library(stringr)
library(dplyr)
library(foreach)
library(stringdist)
library(igraph)
library(data.table)
library(Matrix)
library(doMC)
library(ggrepel)

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
  total_cells_csf   <- length(unique(results[,cellbarcode_col])) 
  total_cells_pb    <- length(unique(results[,cellbarcode_col])) 
  
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
    
    # Calulate percentage of total clones in the patient
    clone_pct     <- (clone_cell_count / total_cells) * 100
    clone_pct_csf <- (clone_csf_count / total_cells_csf) * 100
    clone_pct_pb  <- (clone_pb_count / total_cells_pb) * 100
    
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




### TCR mAb table to GML file for Cytoscape
createTCRnetworkGML <- function(
  data_tcr,
  result_dir      = "./Results",
  output_name     = "cleaned_clusters_all",
  count_cutoff    = as.numeric(2), # UMI count cut-off per contig
  distance_cutoff = as.numeric(0),
  method          = "hamming",
  cores           = 12
){
  registerDoMC(cores=cores)
  print(Sys.time())
  ptm <- proc.time()
  
  # Important files for storing output
  # edge_list_from_distance.R OUTPUT
  overall_cluster_edge_list_file <- file.path(result_dir, paste0("network_edge_list.txt"))
  # write_all_clusters.R OUTPUT
  all_clusters_gml_file <- file.path(result_dir, paste0(output_name, ".gml"))
  all_clusters_csv_file <-file.path(result_dir, paste0(output_name ,".csv"))
  
  ### 1. FORMAT DATA (get_data_bcr.R) -----
  # create target dir
  if (!file.exists(result_dir)) {dir.create(result_dir, recursive=TRUE)}
  
  # UMI filter
  data_tcr <- data_tcr[data_tcr$umi_count>=count_cutoff,]
  
  # Format dataframe so alpha and Beta is on one line
  data_tcr <- formatTCRmabTable(data_tcr)
  
  ### 2. EDGES (edge_list_from_distance.R) -----
  edge_list_final <-foreach(data_of_cluster= split(data_tcr,data_tcr$cluster_criteria), .combine='rbind') %dopar% {
    cluster_by_distance_tcr(data_of_cluster,method, distance_cutoff)
  }
  
  edge_list_final <- unique(edge_list_final)
  
  # ? sort list?
  write.table(edge_list_final, overall_cluster_edge_list_file ,sep="\t")
  
  
  ### 3. WRITE CLUSTERS (write_all_clusters.R) -----
  
  ### Keep Relavent Columns: Updated this line for single cell
  cols_to_keep <- c("Node_ID","cluster_criteria", "COMPARTMENT","STATUS", 
                    "clone_id","junction_aa","CDR3_AA_alpha","SUBTYPE_ANNOTATION",
                    "cell_count","CSF_count","clonotype_CSF_freq","gliph_status","gliph_group")
  vertex<-data_tcr[,cols_to_keep]
  
  vertex <- unique(vertex)
  #names(vertex) <- vertex$Node_ID

  graph_all<- graph.data.frame(edge_list_final, directed=F,vertices=vertex)
  
  write.graph(graph_all, all_clusters_gml_file , format = "gml")
  write.csv(vertex, all_clusters_csv_file )
  
  
  paste(commandArgs(), collapse = " ")
  print(proc.time() - ptm)
}



# Create GLIPH2 network
createGLIPH2clonotypeNetworkGML <- function(
  data_tcr,
  result_dir      = "./Results",
  output_name     = "cleaned_clusters_all",
  count_cutoff    = as.numeric(2), # UMI count cut-off per contig
  distance_cutoff = as.numeric(0),
  method          = "hamming",
  gliphStatusColumn   = "gliph_status",
  gliphInteracColumn = "gliph_group",
  cores           = 12
){
  registerDoMC(cores=cores)
  print(Sys.time())
  ptm <- proc.time()
  
  # Important files for storing output
  # edge_list_from_distance.R OUTPUT
  overall_cluster_edge_list_file <- file.path(result_dir, paste0("network_edge_list.txt"))
  # write_all_clusters.R OUTPUT
  all_clusters_gml_file <- file.path(result_dir, paste0(output_name, ".gml"))
  all_clusters_csv_file <-file.path(result_dir, paste0(output_name ,".csv"))
  
  ### 1. FORMAT DATA (get_data_bcr.R) -----
  # create target dir
  if (!file.exists(result_dir)) {dir.create(result_dir, recursive=TRUE)}
  
  # UMI filter
  data_tcr <- data_tcr[data_tcr$umi_count>=count_cutoff,]
  
  # Format dataframe so alpha and Beta is on one line
  data_tcr <- formatTCRmabTable(data_tcr)
  clone_list <- unique(data_tcr$clone_id)
  
  ### ----Add clonotype-Cell type column
  simple_annot <- data_tcr$SUBTYPE_ANNOTATION
  names(simple_annot) <- data_tcr$clone_id
  #data_tcr$clone_simple_annot <- paste(data_tcr$subtype.cell.annot.manual,data_tcr$MERGED_clonotypeID_imm)
  for (clone in clone_list) {
    clone_sub <- simple_annot[names(simple_annot) %in% clone]
    uni_elements <- length(unique(clone_sub))
    if(uni_elements > 1){
      cd4_cd8 <- paste(as.character(unique(clone_sub)),collapse = " & ")
      simple_annot[names(simple_annot) %in% clone] <- cd4_cd8
    }
  }
  data_tcr$simple_annot <- simple_annot
  
  
  # collapse to clonotype resolution rather than cell resolution
  # -- Remove columns that are not unique (cell barcodes)
  #    Node ID now becomes Clonotype ID
  
  cols_to_omit <- c()
  
  for (clone in clone_list) {
    clone_df <- data_tcr[data_tcr$clone_id %in% clone,]
    for (col in names(clone_df)) {
      one_col <- clone_df[,col]
      uni_elements <- length(unique(one_col))
      if(uni_elements > 1){
        cols_to_omit <- c(cols_to_omit,col)
        if (col == "CDR3_AA_alpha"){
          cat(paste0(clone,"\n")) # sanity check 
        }
      }
    }
  }
  cols_to_omit <- unique(cols_to_omit)
  
  # Collapse data frame to one line per clonotype
  data_tcr_fmt <- data_tcr[,! names(data_tcr) %in% cols_to_omit]
  
  data_tcr_fmt   <- unique(data_tcr_fmt)
  clones_to_keep <-  data_tcr_fmt$clone_id[!is.na(data_tcr_fmt$clone_id)]
  data_tcr_fmt   <- data_tcr_fmt[data_tcr_fmt$clone_id %in% clones_to_keep,]
  row.names(data_tcr_fmt) <- data_tcr_fmt$clone_id
  data_tcr_fmt$Node_ID <- data_tcr_fmt$clone_id
  
  ### 2. EDGES (edge_list_from_distance.R) -----
  
  ## a. Subset GLIPH relavent results
  interac_col        <- grep(gliphInteracColumn,names(data_tcr_fmt))
  interac_status_col <- grep(gliphStatusColumn,names(data_tcr_fmt))
  data_tcr_gliph     <- data_tcr_fmt[data_tcr_fmt[,names(data_tcr_fmt) %in% gliphStatusColumn] %in% "1",]
  interac_df <- data_tcr_gliph[,c("Node_ID","CDR3ab_AA",gliphInteracColumn)]
  
  ## b. Create edge list data frame
  edge_list_final <- data.frame(matrix(nrow = 0,ncol = 2))
  for (i in 1:nrow(interac_df)) {
    one_row  <- interac_df[i,]
    interac1    <- one_row$Node_ID
    interac1_ab <- one_row$CDR3ab_AA
    
    # Find list of all interacting nodes
    interac_grp     <- one_row[,names(one_row) %in% gliphInteracColumn]
    interac_grp_vec <- c()
    if(grepl(":",interac_grp)){
      interac_grp_vec <- as.character(unlist(strsplit(interac_grp, ":")))
    }
    interac_grp_vec <- c(interac_grp_vec,interac_grp)
    
    
    # list all pair interactions for the row
    for (grp in interac_grp_vec) {
      nodes_grp  <- interac_df[interac_df$gliph_group %in% grp,]
      
      # Exclude Interaction 1 from this subset
      nodes_grp   <- nodes_grp[! nodes_grp$CDR3ab_AA %in% interac1_ab,]
      interac2    <- unique(nodes_grp$Node_ID)
      
      # Skip if interaction is not present
      if(length(interac2) == 0){next}
      
      entry       <- data.frame(matrix(nrow = length(interac2), ncol = 2))
      entry$X1 <- interac1
      entry$X2 <- interac2
      edge_list_final <- rbind(edge_list_final,entry)
    }
  }
  
  edge_list_final <- unique(edge_list_final)
  
  write.table(edge_list_final, overall_cluster_edge_list_file ,sep="\t")
  
  
  ### 3. WRITE CLUSTERS (write_all_clusters.R) -----
  
  ### Keep Relavent Columns: Updated this line for single cell
  cols_to_keep <- c("Node_ID", "COMPARTMENT","SAMPLE","STATUS","simple_annot", 
                    "clone_id",gliphStatusColumn,gliphInteracColumn,"cell_count","CSF_count","clonotype_CSF_freq",
                    "CDR3ab_AA")
  vertex<-data_tcr_fmt[,cols_to_keep]
  graph_all<- graph.data.frame(edge_list_final, directed=F,vertices=vertex)
  
  write.graph(graph_all, all_clusters_gml_file , format = "gml")
  write.csv(vertex, all_clusters_csv_file )
  
  
  paste(commandArgs(), collapse = " ")
  print(proc.time() - ptm)
}







