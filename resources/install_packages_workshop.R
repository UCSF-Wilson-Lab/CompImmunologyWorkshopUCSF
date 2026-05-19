# BMS 270 - Immune Repertoire Workshop
#  - Most updated documentation is on github (https://github.com/UCSF-Wilson-Lab/CompImmunologyWorkshopUCSF)

## Set up gfortran and openssl before installing the Bioconductor specific packages

# CRAN R packages
install.packages("Seurat")
install.packages("kableExtra")
install.packages("alakazam")
install.packages("shazam")
install.packages("tigger")
install.packages("dowser")
install.packages("airr")
install.packages("gridExtra")
install.packages("Matrix")
install.packages("ggupset")

install.packages("devtools")
install.packages("knitr")
install.packages("dplyr")
install.packages("stringr")
install.packages("reshape2")
install.packages("gplots")
install.packages("pheatmap")
install.packages("ggplot2")
install.packages("RColorBrewer")
install.packages("stringdist")
install.packages("igraph")
install.packages("data.table")
install.packages("cowplot")
install.packages("ggpubr")
install.packages("enrichplot")
install.packages("ggnewscale")
install.packages("patchwork")


# Bioconductor R packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "devel")

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
BiocManager::install("dreamlet")


# Github repos to install
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

library(devtools)
install_github("UCSF-Wilson-Lab/impactSingleCellToolkit")

# Use older Signac version for installing Azimuth
library(remotes)
remove.packages("Signac")
remotes::install_version("Signac", version = "1.12.0")
remotes::install_github("satijalab/azimuth", ref = "master")
