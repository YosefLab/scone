#' Sum scaling normalization function
#' @details SCONE scaling by library size or summed expression.
#' @export
#' @param ei Numerical matrix. (rows = genes, cols = samples).
#' @return Sum-scaled normalized matrix.
#'
#' @examples
#' ei <- matrix(0:20,nrow = 7)
#' eo <- SUM_FN(ei)
#'
SUM_FN = function (ei) {
  scales = colSums(ei)
  eo = t(t(ei) * mean(scales) / scales)
  return(eo)
}

#' Weighted trimmed mean of M-values (TMM) scaling normalization wrapper
#' function
#' @importFrom edgeR calcNormFactors
#' @details SCONE scaling wrapper for \code{\link[edgeR]{calcNormFactors}}).
#' @export
#' @param ei Numerical matrix. (rows = genes, cols = samples).
#' @return TMM normalized matrix.
#'
#' @examples
#' ei <- matrix(0:20,nrow = 7)
#' eo <- TMM_FN(ei)
#'
TMM_FN = function(ei) {
  size_fac = calcNormFactors(ei, method = "TMM")
  scales = (colSums(ei) * size_fac)
  eo = t(t(ei) * mean(scales) / scales)
  return(eo)
}

#' Relative log-expression (RLE; DESeq) scaling normalization wrapper
#' function
#' @importFrom edgeR calcNormFactors
#' @details SCONE scaling wrapper for \code{\link[edgeR]{calcNormFactors}}).
#' @export
#' @param ei Numerical matrix. (rows = genes, cols = samples).
#' @return RLE normalized matrix.
#'
#' @examples
#' ei <- matrix(0:20,nrow = 7)
#' eo <- DESEQ_FN(ei)
#'
DESEQ_FN = function(ei) {
  size_fac = calcNormFactors(ei, method = "RLE")
  scales = (colSums(ei) * size_fac)
  eo = t(t(ei) * mean(scales) / scales)
  return(eo)
}

#' Upper-quartile (UQ) scaling normalization wrapper function
#' @importFrom edgeR calcNormFactors
#' @details SCONE scaling wrapper for \code{\link[edgeR]{calcNormFactors}}).
#' @export
#' @param ei Numerical matrix. (rows = genes, cols = samples).
#' @return UQ normalized matrix.
#'
#' @examples
#' ei <- matrix(0:20,nrow = 7)
#' eo <- UQ_FN(ei)
#'
UQ_FN = function(ei) {
  size_fac = calcNormFactors(ei, method = "upperquartile")
  scales = (colSums(ei) * size_fac)
  eo = t(t(ei) * mean(scales) / scales)
  return(eo)
}

#' Full-quantile normalization wrapper function
#' @importFrom aroma.light normalizeQuantileRank.matrix
#' @details SCONE "scaling" wrapper for
#'   \code{\link[aroma.light]{normalizeQuantileRank.matrix}}).
#' @export
#' @param ei Numerical matrix. (rows = genes, cols = samples).
#' @return Full-quantile normalized matrix.
#'
#' @examples
#' ei <- matrix(0:20,nrow = 7)
#' eo <- FQ_FN(ei)
#'
FQ_FN = function(ei) {
  eo = normalizeQuantileRank.matrix(ei)
  return(eo)
}

#' @rdname FQ_FN
#' @details Unlike FQ_FN, FQT_FN handles ties carefully (see
#'   \code{\link[limma]{normalizeQuantiles}} for details).
#' @export
#'
#' @examples
#' ei <- matrix(0:20,nrow = 7)
#' eo <- FQT_FN(ei)
#'
FQT_FN = function(ei) {
  eo = normalizeQuantileRank.matrix(ei, ties = TRUE)
  return(eo)
}

#' Centered log-ratio (CLR) normalization wrapper function
#' @importFrom compositions clr
#' @importFrom matrixStats colMedians
#' @details SCONE scaling wrapper for
#'   \code{\link[compositions]{clr}}).
#' @export
#' @param ei Numerical matrix. (rows = genes, cols = samples).
#' @return CLR normalized matrix.
#'
#' @examples
#' ei <- matrix(0:20,nrow = 7)
#' eo <- CLR_FN(ei)
#'
CLR_FN = function (ei)
{
  scale_mat <- t(clr(t(ei))) - log(ei)
  scale_mat[ei == 0] = NA
  scales = exp(-colMedians(scale_mat, na.rm = TRUE))
  eo = t(t(ei) * mean(scales) / scales)
  return(eo)
}

#' Simple deconvolution normalization wrapper
#' @details SCONE scaling wrapper for \code{\link[scran]{computeSumFactors}}).
#' @export
#' @param ei Numerical matrix. (rows = genes, cols = samples).
#' @return scran normalized matrix.
#'
#' @examples
#' ei <- matrix(0:76,nrow = 7)
#' eo <- SCRAN_FN(ei)
#'
SCRAN_FN = function(ei){
  if (!requireNamespace("scran", quietly = TRUE)) {
    stop("scran package needed for SCRAN_FN()")
  }

  scales = scran::calculateSumFactors(ei, sizes = ceiling(sqrt(ncol(ei))))
  eo = t(t(ei) * mean(scales) / scales)
  return(eo)
}

#' PsiNorm normalization wrapper
#' @details SCONE scaling wrapper for \code{\link{PsiNorm}}).
#' @export
#' @param ei Numerical matrix. (rows = genes, cols = samples).
#' @return PsiNorm normalized matrix.
#'
#' @examples
#' ei <- matrix(c(1,0,2,0,2,9,3,0), ncol=2)
#' eo <- PSINORM_FN(ei)
#'
PSINORM_FN = function(ei){
  inv_sf <- pareto.MLE(ei+1)
  eo = t(t(ei) / mean(inv_sf) * inv_sf)
  return(eo)
}
