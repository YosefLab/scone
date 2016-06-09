## SCONE ##

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

### Single-Cell Overview of Normalized Expression data ###

scone is a package to compare and rank the performance of different normalization schemes in real single-cell RNA-seq datasets.

### To Install ###

The easiest way to install the package is to type in R

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("devtools")
biocLite("YosefLab/scone", dependencies=TRUE)
```

