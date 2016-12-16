## SCONE ##

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/YosefLab/scone.svg?branch=master)](https://travis-ci.org/YosefLab/scone)
[![Coverage](https://codecov.io/gh/YosefLab/scone/branch/master/graph/badge.svg)](https://codecov.io/gh/YosefLab/scone)
### Single-Cell Overview of Normalized Expression data ###

SCONE (Single-Cell Overview of Normalized Expression), a package for single-cell RNA-seq data quality control (QC) and normalization. This data-driven framework uses summaries of expression data to assess the efficacy of normalization workflows.

### To Install ###

The easiest way to install the package is to type in R

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("devtools")
biocLite("YosefLab/scone", dependencies=TRUE)
```

Note that SCONE is currently under consideration in Bioconductor and hence requires R-devel (>= 3.4) and Bioconductor devel. 

You can download the latest release of SCONE for R 3.3 [here](https://github.com/YosefLab/scone/releases/tag/v0.99.0).
