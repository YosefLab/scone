## SCONE ##

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/YosefLab/scone.svg?branch=master)](https://travis-ci.org/YosefLab/scone)

### Single-Cell Overview of Normalized Expression data ###

SCONE (Single-Cell Overview of Normalized Expression), a package for single-cell RNA-seq data quality control (QC) and normalization. This data-driven framework uses summaries of expression data to assess the efficacy of normalization workflows.

### To Install ###

The easiest way to install the package is to type in R

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("devtools")
biocLite("YosefLab/scone", dependencies=TRUE)
```

Note that SCONE requires R (>= 3.3) and Bioconductor (>= 3.3). You can download the latest version of R from [CRAN](https://cran.r-project.org). To update Bioconductor to its latest release, follow the instructions [here](https://www.bioconductor.org/install/#update-bioconductor-packages).
