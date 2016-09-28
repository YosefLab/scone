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


## all these should return NULL if the object hasn't have that info
setGeneric(
  name = "get_negconruv",
  def = function(x) {
    standardGeneric("get_negconruv")
  }
)

setGeneric(
  name = "get_negconeval",
  def = function(x) {
    standardGeneric("get_negconeval")
  }
)

setGeneric(
  name = "get_poscon",
  def = function(x) {
    standardGeneric("get_poscon")
  }
)

setGeneric(
  name = "get_qc",
  def = function(x) {
    standardGeneric("get_qc")
  }
)

setGeneric(
  name = "get_bio",
  def = function(x) {
    standardGeneric("get_bio")
  }
)

setGeneric(
  name = "get_batch",
  def = function(x) {
    standardGeneric("get_batch")
  }
)
