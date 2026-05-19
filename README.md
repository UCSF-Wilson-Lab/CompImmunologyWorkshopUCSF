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

This workshop uses **R v4.6.0** and packages installed with **Bioconductor 3.24**. Install instructions for all required R packages and software is provided on Box. In addition, Openssl and GFortran need to be installed. If the below instructions don't work, Hombrew (**brew install**) can also be used to to install Openssl. This option takes extra steps in terms of configuring paths. 

**Install GFortran by downloading the file for your relevant OS**

* Website for all installs [here](https://fortran-lang.org/learn/os_setup/install_gfortran/)
* Mac OS install for **gfortran 14.2** [here](https://github.com/fxcoudert/gfortran-for-macOS/releases)
* If downloading and installing is not working well, you can install with homebrew instead
```{bash}
brew install gcc
```

**Install Openssl in command line for dependencies such as alabastar.base** 
```{bash}
sudo R
source("https://mac.R-project.org/bin/install.R")
install.libs("openssl")
q()
```

```{bash}
# If dependencies are still failing to install because of openssl, then use homebrew
brew install gsl openssl@3
```

* if openssl and or gfortran paths are not being detected while installing R package dependencies, try configuring the paths (Mac OS solution)
```{bash}
# Only run this if other solutions are not working
mkdir -p ~/.R

# Capture the actual OpenSSL and gfortran paths on your machine
OPENSSL_PREFIX="$(brew --prefix openssl@3)"
GFORTRAN_BIN="$(command -v gfortran || true)"

if [ -z "$GFORTRAN_BIN" ]; then
  echo "gfortran is not on PATH."
  echo "Install the current GNU Fortran package for macOS from the R tools page, then rerun."
  exit 1
fi

GFORTRAN_LIBDIR="$(dirname "$("$GFORTRAN_BIN" -print-file-name=libgfortran.dylib)")"

# Makevars config
cat > ~/.R/Makevars <<EOF
OPENSSL_PREFIX := ${OPENSSL_PREFIX}
PKG_CPPFLAGS += -I\$(OPENSSL_PREFIX)/include
PKG_LIBS += -L\$(OPENSSL_PREFIX)/lib

FC := ${GFORTRAN_BIN}
F77 := ${GFORTRAN_BIN}
FLIBS := -L${GFORTRAN_LIBDIR} -lgfortran -lquadmath -lemutls_w -lheapt_w
EOF

# Remove any bad flags from the current shell session
unset LDFLAGS FLIBS CPPFLAGS PKG_LIBS
```


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