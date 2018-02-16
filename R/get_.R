#' @rdname get_params
#'   
#' @param x an object of class \code{\link{SconeExperiment}}.
#'   
#' @return A data.frame containing workflow parameters for each scone workflow.
#'   
#' @export
#' 
setMethod(
  f = "get_params",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    return(x@scone_params)
  })

#' @rdname get_scores
#'   
#' @param x an object of class \code{\link{SconeExperiment}}.
#'   
#' @return \code{get_scores} returns a matrix with all (non-missing) scone 
#'   scores, ordered by average score rank.
#'   
#' @export
#' 
setMethod(
  f = "get_scores",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    scores <- t(na.omit(t(x@scone_scores[,-NCOL(x@scone_scores)])))
    return(scores)
  })

#' @rdname get_scores
#'   
#' @return \code{get_score_ranks} returns a vector of average score ranks.
#'   
#' @export
#' 
setMethod(
  f = "get_score_ranks",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    return(x@scone_scores[,NCOL(x@scone_scores)])
  })


#' @rdname get_negconruv
#'   
#' @param x an object of class \code{\link{SconeExperiment}}.
#'   
#' @return NULL or a logical vector.
#'   
#' @return For \code{get_negconruv} the returned vector indicates which genes
#'   are negative controls to be used for RUV.
#'   
#' @export
#' 
setMethod(
  f = "get_negconruv",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    if(length(x@which_negconruv) == 0) {
      return(NULL)
    } else {
      return(rowData(x)[,x@which_negconruv])
    }
  }
)

#' @rdname get_negconruv
#'   
#' @return For \code{get_negconeval} the returned vector indicates which genes
#'   are negative controls to be used for evaluation.
#'
#' @export
#' 
setMethod(
  f = "get_negconeval",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    if(length(x@which_negconeval) == 0) {
      return(NULL)
    } else {
      return(rowData(x)[,x@which_negconeval])
    }
  }
)

#' @rdname get_negconruv
#'   
#' @return For \code{get_poscon} the returned vector indicates which genes are 
#'   positive controls to be used for evaluation.
#'
#' @export
#'
setMethod(
  f = "get_poscon",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    if(length(x@which_poscon) == 0) {
      return(NULL)
    } else {
      return(rowData(x)[,x@which_poscon])
    }
  }
)

#' @rdname get_qc
#'   
#' @param x an object of class \code{\link{SconeExperiment}}.
#'   
#' @return NULL or the quality control (QC) metric matrix.
#'
#' @export
#'
setMethod(
  f = "get_qc",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    if(length(x@which_qc) == 0) {
      return(NULL)
    } else {
      retval <- as.matrix(colData(x)[, x@which_qc, drop=FALSE])
      return(retval)
    }
  }
)

#' @rdname get_bio
#'   
#' @param x an object of class \code{\link{SconeExperiment}}.
#'   
#' @return NULL or a factor containing bio or batch covariate.
#'
#' @export
#'
setMethod(
  f = "get_bio",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    if(length(x@which_bio) == 0) {
      return(NULL)
    } else {
      return(colData(x)[,x@which_bio])
    }
  }
)

#' @rdname get_bio
#'
#' @export
#'
setMethod(
  f = "get_batch",
  signature = signature(x = "SconeExperiment"),
  definition = function(x) {
    if(length(x@which_batch) == 0) {
      return(NULL)
    } else {
      return(colData(x)[,x@which_batch])
    }
  }
)