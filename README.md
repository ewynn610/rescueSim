# Rescue: REpeated measures Single Cell RNA-seqUEncing Data Simulation
## Getting started 

Rescue is a package for simulating repeated measures scRNA-Seq data. It can be installed from github using the following code:

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
