#' @rdname scone
setGeneric(
  name = "scone",
  def = function(x, ...) {
    standardGeneric("scone")
  }
)

setGeneric(
  name = "get_normalized",
  def = function(x, method) {
    standardGeneric("get_normalized")
  }
)

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
