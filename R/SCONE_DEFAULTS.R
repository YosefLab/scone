#' Upper-quartile normalization wrapper.
#' @importFrom EDASeq betweenLaneNormalization
#' @export
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return Upper-quartile normalized matrix (scaled by sample UQ)
UQ_FN = function(ei){
  eo = betweenLaneNormalization(ei, which="upper", round = FALSE)
  return(eo)
}

uq_f <- function(which="upper", round = FALSE){
  function(x) {
    betweenLaneNormalization(x, which = which, round = round)
  }
}

# setMethod(f="UQ_FN",
#          signature="SummarizedExperiment",
#          definition=function(x) {
#            exprs(x) <- UQ_FN(exprs(x))
#            return(x)
#          })
#
# setMethod(f="UQ_FN",
#           signature="matrix",
#           definition= function(ei){
#             eo = betweenLaneNormalization(ei, which="upper", round = FALSE)
#             return(eo)
#           })

#' Upper-quartile normalization applied to positive data.
#' @export
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return Upper-quartile (positive) normalized matrix (scaled by sample UQ)
UQ_FN_POS = function(ei){
  zei = ei
  zei[ei == 0] = NA
  q = apply(zei, 2, quantile, 0.75, na.rm = TRUE)
  zeo = t(t(zei)/q)*mean(q)
  eo = zeo
  eo[ei == 0] = 0
  return(eo)
}

#' Full-Quantile normalization wrapper.
#' @importFrom aroma.light normalizeQuantileRank.matrix
#' @export
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return FQ normalized matrix.
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

#' Full-Quantile normalization applied to positive data.
#' @export
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return FQ (positive) normalized matrix.
FQ_FN_POS = function(ei){

  # Vector of integers used for computation
  base_rank = 1:nrow(ei)

  # Quantile Index Matrix: Values between 0 and 1 corresponding to quantile
  quant_mat = NULL
  # Re-ordered Data Matrix
  x_mat = NULL
  # For each sample:
  for (i in 1:dim(ei)[2]){
    # Sort data and replace zeroes with NA
    x_mat = cbind(x_mat,rev(sort(ei[,i])))
    x_mat[x_mat == 0 ] = NA
    # Compute corresponding quantile indices for that sample
    quant = base_rank/sum(ei[,i]>0)
    quant[quant > 1] = NA
    quant_mat = cbind(quant_mat,quant)
  }

  # Vector form of quantile index matrix
  quant_out = as.numeric(as.vector(quant_mat))
  # Interpolation Matrix (Values of all quantiles)
  inter_mat = rep(0,length(quant_out))
  ob_counts = rep(0,length(quant_out)) # Number of observations for averaging
  # For each sample
  for (i in 1:dim(ei)[2]){
    # Produce spline interpolation for that sample, value ~ quantile index
    x1 = na.omit(quant_mat[,i])
    y1 = na.omit(x_mat[,i])
    # Evaluated at all possible quantile indices
    inter = approx(x1,y1,xout = quant_out, rule = 2)$y
    ob_counts = ob_counts + !is.na(inter)
    inter[is.na(inter)] = 0
    inter_mat = inter_mat + inter
  }

  # Average over the interpolated values from all samples
  inter_mean = inter_mat/ob_counts

  ## Substituting Mean Interpolated Values for Expression Values and Return
  eo = matrix(inter_mean,ncol = dim(ei)[2])
  eo[is.na(eo)] = 0
  for (i in 1:dim(ei)[2]){
    eo[,i] = rev(eo[,i])[order(order(ei[,i]))]
  }

  rownames(eo) = rownames(ei)
  colnames(eo) = colnames(ei)
  return(eo)
}

#' DESeq normalization wrapper.
#' @importFrom DESeq estimateSizeFactorsForMatrix
#' @export
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return DESeq normalized matrix (scaled by sample )
DESEQ_FN = function(ei){
  size_fac = estimateSizeFactorsForMatrix(ei)
  eo = t(t(ei)/size_fac)
  return(eo)
}

#' DESeq normalization applied to positive data.
#' @export
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return DESeq (positive) normalized matrix (scaled by sample )
DESEQ_FN_POS = function(ei){
  if(any(ei < 0)){
    stop("Negative values in input.")
  }
  if(!is.null(dim(ei))){
    y = ei
    y[y == 0] = NA # Matrix with zeroes replaced w/ NA
    geom_mean = exp(apply(log(y),1,sum,na.rm = TRUE)/rowSums(!is.na(y))) # Compute Geometric Mean of Expression for Each Gene (Use positive data only)
  }else{
    stop("Null imput matrix dimension.")
  }
  if(!any(geom_mean > 0)){
    stop("Geometric mean non-positive for all genes.")
  }
  ratios = ei / geom_mean # Divide each Expression Value by Geometric Mean of Corresponding Gene
  if(any(is.infinite(geom_mean))){
    stop("Infinite mean! This should never happen :-<")
  }
  ratios = ratios[geom_mean > 0,] # Ignore Genes with Zero Mean
  y = ratios
  y[y == 0] = NA
  size = apply(y,2,median,na.rm = TRUE) # Size factor taken as median of ratios (positive data only)
  if(any(size == 0)){
    stop("Zero library size. This should never happen :-(")
  }
  eo = t(t(ei)/size)*mean(size)
  return(eo)
}

#' TMM normalization wrapper.
#' @importFrom edgeR calcNormFactors
#' @export
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return TMM normalized matrix (scaled by sample )
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
