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
#'   object adding an assay with the normalized count matrix. If the input
#'   is a matrix-like object PsiNorm returns a matrix with the
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
#' @importFrom DelayedMatrixStats colSums2 colMins
#' @importFrom sparseMatrixStats colSums2 colMins
#' @importClassesFrom SparseArray SparseMatrix
#' @importFrom SparseArray nzwhich
setMethod(
  f = "PsiNorm",
  signature = signature(x = "ANY"),
  definition = function(x){
    ## Uses only two transpostions (instead of four in original implementation).
    ## Preserves sparsity.
    tx <- t(x)
    sf <- compute_transposed_PsiNorm_sf(tx)
    t(mat_v_div(tx, sf))
})

## Preserves sparsity.
compute_transposed_PsiNorm_sf <- function(tx) {
  ## Temporary workaround until operations from the Math group (e.g. log1p)
  ## are supported on a SparseMatrix object of type "integer".
  ## See https://github.com/Bioconductor/SparseArray/blob/devel/TODO
  if (type(tx) != "double")
    type(tx) <- "double"
  x2 <- log1p(tx)
  m2 <- MatrixGenerics::rowMins(x2)
  x2 <- mat_v_sub(x2, m2)
  MatrixGenerics::rowSums2(x2) / ncol(tx)
}

## Preserves sparsity.
computePsiNormSF <- function(x) {
  compute_transposed_PsiNorm_sf(t(x))
}

.build_linear_index <- function(nr, nc, i) {
  stopifnot(is.integer(nr), length(nr) == 1L,
            is.integer(nc), length(nc) == 1L, is.integer(i))
  rep.int(nr * (seq_len(nc) - 1L), rep.int(length(i), nc)) + i
}

## Implements optimized substraction between a matrix-like object 'mat' and
## an ordinary vector 'v' with mostly zeros that accepts a SparseMatrix
## object (in which case it also returns a SparseMatrix object). Note that
## the SparseArray package does not support this operation out-of-the-box at
## the moment.
mat_v_sub <- function(mat, v) {
  mat_dim <- dim(mat)
  stopifnot(length(mat_dim) == 2L, is.vector(v), length(v) == mat_dim[[1L]])
  nzidx <- nzwhich(v)  # indices of nonzero values in 'v'
  if (length(nzidx) == 0L)
    return(mat)  # no-op
  if (!is(mat, "SparseMatrix"))
    return(mat - v)
  ## Works well if 'nzidx' is short (i.e. 'v' contains mostly zeros).
  Lidx <- .build_linear_index(mat_dim[[1L]], mat_dim[[2L]], nzidx)
  mat[Lidx] <- mat[Lidx] - v[nzidx]
  mat
}

## Implements division between a matrix-like object 'mat' and an ordinary
## vector 'v' with mostly nonzero/non-NA values that accepts a SparseMatrix
## object (in which case it also returns a SparseMatrix object). Note that
## the SparseArray package only supports this operation when 'v' contains
## no zeros or NAs at the moment.
mat_v_div <- function(mat, v) {
  mat_dim <- dim(mat)
  stopifnot(length(mat_dim) == 2L, is.vector(v), length(v) == mat_dim[[1L]])
  idx <- which(is.na(v) | v == 0)
  if (length(idx) == 0L || !is(mat, "SparseMatrix"))
    return(mat / v)
  ## Works well if 'idx' is short (i.e. if 'v' contains very few zeros/NAs).
  v2 <- v
  v2[idx] <- 1L
  mat <- mat / v2
  Lidx <- .build_linear_index(mat_dim[[1L]], mat_dim[[2L]], idx)
  mat[Lidx] <- mat[Lidx] / v[idx]
  mat
}

