#' Upper-quartile normalization wrapper.
#' @importFrom EDASeq betweenLaneNormalization
#' @details SCONE scaling wrapper for \code{\link[EDASeq]{betweenLaneNormalization}}).
#' @export
#' @param ei Numerical matrix. (rows = genes, cols = samples).
#' @return Upper-quartile normalized matrix.
UQ_FN = function(ei){
  eo = betweenLaneNormalization(ei, which="upper", round = FALSE)
  return(eo)
}

#' Upper-quartile normalization derived from positive data.
#' @details SCONE scaling function scales expression data by upper quartile of positive data.
#' @export
#' @param ei Numerical matrix. (rows = genes, cols = samples).
#' @return Upper-quartile (positive) normalized matrix.
UQ_FN_POS = function(ei){
  zei = ei
  is_zero = (ei == 0)
  zei[is_zero] = NA
  
  q = apply(zei, 2, quantile, 0.75, na.rm = TRUE)
  zeo = t(t(zei)/q)*mean(q)
  
  eo = zeo
  eo[is_zero] = 0
  return(eo)
}

#' Full-quantile normalization wrapper.
#' @importFrom aroma.light normalizeQuantileRank.matrix
#' @details SCONE "scaling" wrapper for \code{\link[aroma.light]{normalizeQuantileRank.matrix}}).
#' @export
#' @param ei Numerical matrix. (rows = genes, cols = samples).
#' @return Full-quantile normalized matrix.
FQ_FN = function(ei){
  eo = normalizeQuantileRank.matrix(ei)
  return(eo)
}

#' @rdname FQ_FN
#' @details FQT_FN handles ties carefully (see \code{\link[limma]{normalizeQuantiles}}).
#' @export
FQT_FN = function(ei){
  eo = normalizeQuantileRank.matrix(ei, ties = TRUE)
  return(eo)
}

#' DESeq size factor normalization wrapper.
#' @importFrom DESeq estimateSizeFactorsForMatrix
#' @details SCONE scaling wrapper for \code{\link[DESeq]{estimateSizeFactorsForMatrix}}).
#' @export
#' @param ei Numerical matrix. (rows = genes, cols = samples).
#' @return DESeq size factor normalized matrix.
DESEQ_FN = function(ei){
  size_fac = estimateSizeFactorsForMatrix(ei)
  eo = t(t(ei)/size_fac)
  return(eo)
}

#' DESeq size factor normalization derived from positive data.
#' @details SCONE scaling function scales expression data by DESeq size factor derived from positive data.
#' @export
#' @param ei Numerical matrix. (rows = genes, cols = samples).
#' @return DESeq size factor (positive) normalized matrix.
DESEQ_FN_POS = function(ei){
  
  if(any(ei < 0)){stop("Negative values in input.")}
  
  if(!is.null(dim(ei))){
    y = ei
    y[y == 0] = NA # Matrix with zeroes replaced w/ NA
    geom_mean = exp(apply(log(y),1,sum,na.rm = TRUE)/rowSums(!is.na(y))) # Compute Geometric Mean of Expression for Each Gene (Use positive data only)
  }else{stop("Null imput matrix dimension.")}
  if(!any(geom_mean > 0)){stop("Geometric mean non-positive for all genes.")}
  
  # Divide each Expression Value by Geometric Mean of Corresponding Gene
  ratios = ei / geom_mean
  if(any(is.infinite(geom_mean))){stop("Infinite mean for some gene.")}
  
  # Ignore Genes with Zero Geometric Mean
  ratios = ratios[geom_mean > 0,]
  
  # Size factor taken as median of ratios (positive data only)
  y = ratios
  y[y == 0] = NA
  size = apply(y,2,median,na.rm = TRUE) 
  if(any(size == 0)){stop("Zero library size for some sample.")}
  
  eo = t(t(ei)/size)*mean(size)
  return(eo)
}

#' Weighted trimmed mean of M-values (TMM) normalization wrapper.
#' @importFrom edgeR calcNormFactors
#' @details SCONE scaling wrapper for \code{\link[edgeR]{calcNormFactors}}).
#' @export
#' @param ei Numerical matrix. (rows = genes, cols = samples).
#' @return TMM normalized matrix.
TMM_FN = function(ei){
  size_fac = calcNormFactors(ei,method = "TMM")
  eo = t(t(ei)/size_fac)
  return(eo)
}

#' LSF normalization wrapper.
#' @importFrom scran computeSumFactors
#' @importFrom scran quickCluster
#' @export
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return LSF normalized matrix (scaled by sample )
LSF_FN = function(ei){
  clusters <- quickCluster(ei, min.size = 20)
  size_fac = computeSumFactors(ei, cluster=clusters, sf.out = TRUE)
  eo = t(t(ei)/size_fac)
  return(eo)
}
