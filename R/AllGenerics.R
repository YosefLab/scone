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
#' @param x a \code{\link{sconeExperiment}} object containing the results of
#'   \code{\link{scone}}.
#' @param method character or numeric. Either a string identifying the
#'   normalization scheme to be retrieved, or a numeric index with the rank of
#'   the normalization method to retrieve (according to scone ranking of
#'   normalizations).
#' @param ... additional arguments for specific methods.
#'
#' @return A matrix of normalized counts in log-scale.
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
#' @param x a \code{\link{sconeExperiment}} object containing the results of
#'   \code{\link{scone}}.
#' @param method character or numeric. Either a string identifying the
#'   normalization scheme to be retrieved, or a numeric index with the rank of
#'   the normalization method to retrieve (according to scone ranking of
#'   normalizations).
#'
#' @return The design matrix.
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
#' @param methods either character or numeric specifying the normalizations to select.
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
#' @aliases get_qc,SconeExperiment-method
setGeneric(
  name = "get_qc",
  def = function(x) {
    standardGeneric("get_qc")
  }
)

#' Get Factor of Biological Conditions and Batch
#'
#' @aliases get_batch get_bio,SconeExperiment-method
#'   get_batch,SconeExperiment-method
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
#' @aliases get_scores get_score,SconeExperiment-method get_score_ranks
#'   get_score_ranks,SconeExperiment-method
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
setGeneric(
  name = "get_params",
  def = function(x) {
    standardGeneric("get_params")
  }
)
