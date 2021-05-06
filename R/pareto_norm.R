#' PsiNorm: scaling normalization based on the Pareto distribution
#'
#' Normalization of a raw counts matrix using the estimate of the shape
#' parameter of the Pareto distribution.
#'
#' @param x A SingleCellExperiment/SummarizedExperiment object or a
#'   matrix=like object with genes in rows and samples in columns.
#' @param whichAssay if x is a SingleCellExperiment/SummarizedExperiment the
#'   assay with the counts to normalize (default to 1).
#' @param assayName if x is a SummarizedExperiment the name of the assay in
#'   which to save the normalized data (default to "PsiNorm").
#' @importFrom MatrixGenerics colMins colSums2
#' @import SummarizedExperiment
#' @importFrom methods setMethod
#'
#' @return If the input is a SingleCellExperiment object the function returns
#'   the same object adding as sizeFactors those computed by PsiNorm. If the
#'   object is a SummarizedExperiment object, the function returns the same
#'   objject adding an assay with the normalized count matrix. If the input
#'   is a matrix-like object pareto_norm returns a matrix with the
#'   same dimensions containing the normalized counts.
#' @export
#' @author Matteo Borella and Davide Risso
#'
#' @rdname PsiNorm
#' @examples
#' m<-matrix(c(1,0,2,0,2,9,3,0), ncol=2)
#' sce<-SingleCellExperiment::SingleCellExperiment(assays=list(counts=m))
#'
#' sce<-PsiNorm(sce) # SingleCellExperiment object
#' norm.matrix<-PsiNorm(m) # normalized matrix object
#'
setMethod(
  f = "PsiNorm",
  signature = signature(x = "SummarizedExperiment"),
  definition = function(x, whichAssay = 1, assayName = "PsiNorm"){
    assay(x, assayName) <- PsiNorm(assay(x, whichAssay))
    return(x)
})

#' @rdname PsiNorm
#' @export
#' @import SingleCellExperiment
setMethod(
  f = "PsiNorm",
  signature = signature(x = "SingleCellExperiment"),
  definition = function(x, whichAssay = "counts"){
    sf <- computePsiNormSF(assay(x, whichAssay))
    sizeFactors(x) <- sf
    return(x)
})

#' @rdname PsiNorm
#' @export
setMethod(
  f = "PsiNorm",
  signature = signature(x = "ANY"),
  definition = function(x){
    sf <- computePsiNormSF(x)
    t(t(x)/sf)
})

# can we simplify the function and assume min is always 1?
pareto.MLE <- function(sce) {
  n <- nrow(sce)
  m <- MatrixGenerics::colMins(sce)
  a <- n/MatrixGenerics::colSums2(t(t(log(sce)) - log(m)))
  return(a)
}

computePsiNormSF <- function(x) {
  1/pareto.MLE(x+1)
}
