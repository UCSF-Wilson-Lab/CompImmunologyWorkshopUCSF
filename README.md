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
  │   ├── 4_find_candidates.Rmd
  │   └── 5_vibe_vdj_networks.Rmd (optional)
  │
  ├── 4_FindCandidates
  │   └── all_CSF_GLIPH2
  ├── data
  │   ├── input_immcantation
  │   ├── results_immcantation
  │   └── input_antigen_discovery
  ├── figures
  ├── objects
  └── resources
      ├── functions_workshop.R
      └── functions_pub_networks.R
```

## Dependencies

This workshop uses **R v4.6.0** and packages installed with **Bioconductor 3.24 or 3.23**. Install instructions for all required R packages and software is provided on Box. 

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

**Install an older version of Signac to avoid errors with Azimuth installation**
```{r}
library(remotes)
remove.packages("Signac")
remotes::install_version("Signac", version = "1.12.0")
remotes::install_github("satijalab/azimuth", ref = "master")
```

**Github repos to install:**

* [Azimuth](https://github.com/satijalab/azimuth "Azimuth")
* [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder "DoubletFinder")
* [impactSingleCellToolkit](https://github.com/UCSF-Wilson-Lab/impactSingleCellToolkit "impactSingleCellToolkit")