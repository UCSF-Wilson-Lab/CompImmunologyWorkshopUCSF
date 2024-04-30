# BMS 270 Computational Immunology Workshop UCSF

This codebase covers an analysis of single cell RNA-Seq and VDJ data as part of the course BMS 270. All input data and results are provided through Box prior to the workshop.

## Structure

* Both single cell RNA-Seq and VDJ data was pre-processed ahead of time and run through 10x's cellranger software. Afterwards, filtered contigs outputted from cellranger were then inputted into Immcantation and the output was stored in the **data** sub-directory. 
* All data was combined and pre-processed (```0_preprocessing.Rmd```) to omit doublets and annotate cell types. These processed data were stored in the **objects** directory and serves as the input for the workshop. 
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



