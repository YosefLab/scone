## SCONE ##

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/scone.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/scone)
[![R-CMD-check](https://github.com/YosefLab/scone/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/YosefLab/scone/actions)
<!-- badges: end -->

### Single-Cell Overview of Normalized Expression data ###

SCONE (Single-Cell Overview of Normalized Expression), a package for single-cell RNA-seq data quality control (QC) and normalization. This data-driven framework uses summaries of expression data to assess the efficacy of normalization workflows.

### Install from Bioconductor ###

We recommend installation of the package via bioconductor.

```{r}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("scone")
```

### Install from Github ###

Usually not recommended. To download the development version of the package, use

```{r}
BiocManager::install("YosefLab/scone")
```

### Install for R 3.3 ###

You can download the latest release of SCONE for R 3.3 [here](https://github.com/YosefLab/scone/releases/tag/v0.99.0).
This is useful only for reproducing old results.
