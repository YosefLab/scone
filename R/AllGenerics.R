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

setGeneric(
  name = "select_methods",
  def = function(x, index) {
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
