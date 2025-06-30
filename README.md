# rescueSim: Repeated measures single cell RNA-sequencing data simulation

**rescueSim** is an R package for simulating repeated measures (paired/longitudinal) single-cell RNA-sequencing (scRNA-seq) data.  
It is designed to support power analysis, benchmarking, and method development for studies involving complex scRNA-seq designs.

## Getting started 

### Install dependencies

```{r}
## CRAN packages
install.packages(c("checkmate", "dplyr", "gtools", "Matrix"))

## Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("scater", "edgeR", "scran", "scater", "scuttle", "SingleCellExperiment", "MAST"))
```
### Install rescueSim from GitHub
**Note:** The package vignettes are available on the package website. If you'd prefer to build the vignettes locally during installation, you can set `build_vignettes = TRUE`, though this may take longer.

```{r}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("https://github.com/ewynn610/rescueSim",  build_vignettes = FALSE)
```
