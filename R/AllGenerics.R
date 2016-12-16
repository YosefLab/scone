#' @rdname scone
setGeneric(
  name = "scone",
  def = function(x, ...) {
    standardGeneric("scone")
  }
)

#' Retrieve Normalized Matrix
#' 
#' Given a \code{SconeExperiment} object created by a call to scone, it will 
#' return a matrix of normalized counts (in log scale if \code{log=TRUE}).
#' 
#' @details If \code{\link{scone}} was run with \code{return_norm="in_memory"},
#'   this function simply retrieves the normalized data from the \code{assays} 
#'   slote of \code{object}.
#'   
#' @details If \code{\link{scone}} was run with \code{return_norm="hdf5"}, this
#'   function will read the normalized matrix from the specified hdf5 file.
#'   
#' @details If \code{\link{scone}} was run with \code{return_norm="no"}, this 
#'   function will compute the normalized matrix on the fly.
#'   
#' @param x a \code{\link{SconeExperiment}} object containing the results of 
#'   \code{\link{scone}}.
#' @param method character or numeric. Either a string identifying the 
#'   normalization scheme to be retrieved, or a numeric index with the rank of 
#'   the normalization method to retrieve (according to scone ranking of 
#'   normalizations).
#' @param ... additional arguments for specific methods.
#'   
#' @return A matrix of normalized counts in log-scale.
#'   
#' @examples
#' set.seed(42)
#' mat <- matrix(rpois(500, lambda = 5), ncol=10)
#' colnames(mat) <- paste("X", 1:ncol(mat), sep="")
#' obj <- SconeExperiment(mat)
#' res <- scone(obj, scaling=list(none=identity, uq=UQ_FN),
#'            evaluate=TRUE, k_ruv=0, k_qc=0, 
#'            eval_kclust=2, bpparam = BiocParallel::SerialParam())
#' top_norm = get_normalized(res,1)
#'            
#' 
setGeneric(
  name = "get_normalized",
  def = function(x, method, ...) {
    standardGeneric("get_normalized")
  }
)

#' Retrieve Design Matrix
#'
#' Given a \code{SconeExperiment} object created by a call to scone, it will
#' return the design matrix of the selected method.
#'
#' @param x a \code{\link{SconeExperiment}} object containing the results of
#'   \code{\link{scone}}.
#' @param method character or numeric. Either a string identifying the
#'   normalization scheme to be retrieved, or a numeric index with the rank of
#'   the normalization method to retrieve (according to scone ranking of
#'   normalizations).
#'
#' @return The design matrix.
#' 
#' @examples
#' set.seed(42)
#' mat <- matrix(rpois(500, lambda = 5), ncol=10)
#' colnames(mat) <- paste("X", 1:ncol(mat), sep="")
#' obj <- SconeExperiment(mat, bio = factor(rep(c(1,2),each = 5)),
#'            batch = factor(rep(c(1,2),times = 5)))
#' res <- scone(obj, scaling=list(none=identity, uq=UQ_FN),
#'            evaluate=TRUE, k_ruv=0, k_qc=0, 
#'            adjust_batch = "yes", adjust_bio = "yes",
#'            eval_kclust=2, bpparam = BiocParallel::SerialParam())
#' design_top = get_design(res,1)
#' 
setGeneric(
  name = "get_design",
  def = function(x, method) {
    standardGeneric("get_design")
  }
)

#' Get a subset of normalizations from a SconeExperiment object
#' 
#' @description This method let a user extract a subset of normalizations. This
#'   is useful when the original dataset is large and/or many normalization 
#'   schemes have been applied.
#'   
#' @description In such cases, the user may want to run scone in mode 
#'   \code{return_norm = "no"}, explore the results, and then select the top 
#'   performing methods for additional exploration.
#'   
#' @param x a \code{SconeExperiment} object.
#' @param methods either character or numeric specifying the normalizations to
#'   select.
#'   
#' @return A \code{SconeExperiment} object with selected method data.
#'   
#' @examples
#' set.seed(42)
#' mat <- matrix(rpois(500, lambda = 5), ncol=10)
#' colnames(mat) <- paste("X", 1:ncol(mat), sep="")
#' obj <- SconeExperiment(mat)
#' res <- scone(obj, scaling=list(none=identity, uq=UQ_FN),
#'            evaluate=TRUE, k_ruv=0, k_qc=0, 
#'            eval_kclust=2, bpparam = BiocParallel::SerialParam())
#' select_res = select_methods(res,1:2)
#' 
setGeneric(
  name = "select_methods",
  def = function(x, methods) {
    standardGeneric("select_methods")
  }
)

#' Get Negative and Positive Controls
#' 
#' @aliases get_negconeval get_poscon get_negconruv,SconeExperiment-method 
#'   get_negconeval,SconeExperiment-method get_poscon,SconeExperiment-method
#' 
#' @examples
#' set.seed(42)
#' mat <- matrix(rpois(500, lambda = 5), ncol=10)
#' colnames(mat) <- paste("X", 1:ncol(mat), sep="")
#' obj <- SconeExperiment(mat,negcon_ruv = 1:50 %in% 1:10,
#'            negcon_eval = 1:50 %in% 11:20,
#'            poscon = 1:50 %in% 21:30)
#' negcon_ruv = get_negconruv(obj)
#' negcon_eval = get_negconeval(obj)
#' poscon = get_poscon(obj)
#'
setGeneric(
  name = "get_negconruv",
  def = function(x) {
    standardGeneric("get_negconruv")
  }
)

#' @rdname get_negconruv
setGeneric(
  name = "get_negconeval",
  def = function(x) {
    standardGeneric("get_negconeval")
  }
)

#' @rdname get_negconruv
setGeneric(
  name = "get_poscon",
  def = function(x) {
    standardGeneric("get_poscon")
  }
)

#' Get Quality Control Matrix
#' 
#' @examples
#' set.seed(42)
#' mat <- matrix(rpois(500, lambda = 5), ncol=10)
#' colnames(mat) <- paste("X", 1:ncol(mat), sep="")
#' obj <- SconeExperiment(mat,
#'          qc = cbind(colSums(mat),colSums(mat > 0)))
#' qc = get_qc(obj)
#'
#' @aliases get_qc,SconeExperiment-method
setGeneric(
  name = "get_qc",
  def = function(x) {
    standardGeneric("get_qc")
  }
)

#' Get Factor of Biological Conditions and Batch
#'
#' @aliases get_bio get_batch get_bio,SconeExperiment-method
#'   get_batch,SconeExperiment-method
#'   
#' @examples 
#' set.seed(42)
#' mat <- matrix(rpois(500, lambda = 5), ncol=10)
#' colnames(mat) <- paste("X", 1:ncol(mat), sep="")
#' obj <- SconeExperiment(mat, bio = factor(rep(c(1,2),each = 5)),
#'            batch = factor(rep(c(1,2),times = 5)))
#' bio = get_bio(obj)
#' batch = get_batch(obj)
#' 
setGeneric(
  name = "get_bio",
  def = function(x) {
    standardGeneric("get_bio")
  }
)

#' @rdname get_bio
setGeneric(
  name = "get_batch",
  def = function(x) {
    standardGeneric("get_batch")
  }
)

#' Extract scone scores
#' 
#' @aliases get_scores get_scores,SconeExperiment-method get_score_ranks 
#'   get_score_ranks,SconeExperiment-method
#'   
#' @examples
#' set.seed(42)
#' mat <- matrix(rpois(500, lambda = 5), ncol=10)
#' colnames(mat) <- paste("X", 1:ncol(mat), sep="")
#' obj <- SconeExperiment(mat)
#' res <- scone(obj, scaling=list(none=identity, uq=UQ_FN),
#'            evaluate=TRUE, k_ruv=0, k_qc=0, 
#'            eval_kclust=2, bpparam = BiocParallel::SerialParam())
#' scores = get_scores(res)
#' score_ranks = get_score_ranks(res)
#' 
setGeneric(
  name = "get_scores",
  def = function(x) {
    standardGeneric("get_scores")
  }
)

#' @rdname get_scores
setGeneric(
  name = "get_score_ranks",
  def = function(x) {
    standardGeneric("get_score_ranks")
  }
)

#' Extract scone parameters
#' 
#' @aliases get_params get_params,SconeExperiment-method
#' 
#' @examples
#' set.seed(42)
#' mat <- matrix(rpois(500, lambda = 5), ncol=10)
#' colnames(mat) <- paste("X", 1:ncol(mat), sep="")
#' obj <- SconeExperiment(mat)
#' res <- scone(obj, scaling=list(none=identity, uq=UQ_FN),
#'            run = FALSE, k_ruv=0, k_qc=0, eval_kclust=2)
#' params = get_params(res)
#' 
setGeneric(
  name = "get_params",
  def = function(x) {
    standardGeneric("get_params")
  }
)
