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
  │   └── 4_find_candidates.Rmd
  │
  ├── 4_FindCandidates
  │   └── all_CSF_GLIPH2
  ├── data
  │   ├── input_immcantation
  │   └── results_immcantation
  ├── objects
  └── resources
```

## Dependencies

This workshop supports **R v4.3.2** or later. Install instructions for all required R packages and software is provided on Box. 

**CRAN R packages to install:**
```{r}

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

**Github repos to install:**

* [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder "DoubletFinder")
* [impactSingleCellToolkit](https://github.com/UCSF-Wilson-Lab/impactSingleCellToolkit "impactSingleCellToolkit")