# RESCUE: Repeated measures single cell RNA-sequencing data simulation
## Getting started 

RESCUE is a package for simulating repeated measures scRNA-Seq data. To install the package, start by ensuring all dependenciesare installed:
```{r}
## Install CRAN dependencies
install.packages(c("checkmate", "dplyr", "gtools", "Matrix"))

## Install bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("scater", "edgeR", "scran", "scater", "scuttle", "SingleCellExperiment"))

```

The RESCUE package can be installed from github using the following code:
```{r}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("https://github.com/ewynn610/Rescue",  build_vignettes = T)
```

Once the package is installed, to access a vignette introducing the package workflow use:
```{r}
vignette("rescueVignette")
```
