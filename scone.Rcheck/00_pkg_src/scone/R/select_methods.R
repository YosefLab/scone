#' @describeIn select_methods 
#'   If 
#'   \code{methods} is a character, it will return the subset of
#'   methods named in \code{methods} (only perfect match). The 
#'   string must be a subset of the \code{row.names} of the slot
#'   \code{scone_params}.
#'   
#' @export
#' 
setMethod(
  f = "select_methods",
  signature = signature(x = "SconeExperiment", methods = "character"),
  definition =  function(x, methods) {

    if(!all(methods %in% rownames(x@scone_params))) {
      stop("`methods` must be part of the row.names of the slot scone_params.")
    }

    retval <- x

    if(x@scone_run == "in_memory") {
      assays(retval) <- assays(x)[methods]
    }

    retval@scone_params <- x@scone_params[methods,]

    if(!all(is.na(x@scone_scores))) {
      retval@scone_scores <- x@scone_scores[methods,]
      retval@scone_metrics <- x@scone_metrics[methods,]
    }

    validObject(retval)
    return(retval)
  }
)

#' @describeIn select_methods 
#'  If
#'  \code{methods} is a numeric, it will return the subset of methods
#'  according to the scone ranking.
#'   
#' @details The numeric method will always return the normalization 
#'   corresponding to the \code{methods} rows of the \code{scone_params} slot. 
#'   This means that if \code{\link{scone}} was run with \code{eval=TRUE}, 
#'   \code{select_methods(x, 1:3)} will return the top three ranked method. If 
#'   \code{\link{scone}} was run with \code{eval=FALSE}, it will return the 
#'   first three normalization in the order saved by scone.
#'   
#' @export
#' 
setMethod(
  f = "select_methods",
  signature = signature(x = "SconeExperiment", methods = "numeric"),
  definition =  function(x, methods) {

    norm_methods <- rownames(x@scone_params)[methods]

    return(select_methods(x, norm_methods))
  }
)

