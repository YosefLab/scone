#' @describeIn get_design 
#'   If 
#'   \code{method} is a character, it will return the design
#'   matrix corresponding to the normalization scheme specified
#'   by the character string. The string must be one of the
#'   \code{row.names} of the slot \code{scone_params}.
#'   
#' @export
#' 
#' 
setMethod(
  f = "get_design",
  signature = signature(x = "SconeExperiment", method = "character"),
  definition =  function(x, method) {

    params <- unlist(x@scone_params[method,])

    ruv_factors <- qc_factors <- NULL
    if(params[3] != "no_uv") {

      k <- as.numeric(strsplit(params[3], "=")[[1]][2])

      if(grepl("ruv", params[3])) {

        if(length(x@impute_args) > 0) {
          impute_args <- unlist(x@impute_args)
        } else {
          impute_args <- NULL
        }

        imputed <- x@imputation_fn[[params[1]]](assay(x), impute_args)
        scaled <- x@scaling_fn[[params[2]]](imputed)

        if(x@rezero) {
          toz = assay(x) <= 0
          scaled <- scaled - scaled*toz
        }

        r <- RUVg(log1p(scaled), get_negconruv(x), k, isLog=TRUE)$W
        ruv_factors <- list(r)
        names(ruv_factors) <- paste(params[1:2], collapse="_")
      }

      if(grepl("qc", params[3])) {
        qc_factors <- prcomp(get_qc(x), center=TRUE, scale=TRUE)$x
      }

    }

    parsed <- .parse_row(params, get_bio(x), get_batch(x), 
                         ruv_factors, qc_factors)
    design_mat <- make_design(parsed$bio, parsed$batch, parsed$W,
                              nested=(x@nested &
                                        !is.null(parsed$bio) &
                                        !is.null(parsed$batch)))

    return(design_mat)
  }
)

#' @describeIn get_design
#'   If
#'   \code{method} is a numeric, it will return the design matrix
#'   according to the scone ranking.
#'   
#' @details The numeric method will always return the design matrix 
#'   corresponding to row \code{method} of the \code{scone_params} 
#'   slot. This means that if \code{\link{scone}} was run with 
#'   \code{eval=TRUE}, \code{get_design(x, 1)} will return the top
#'   ranked method. If \code{\link{scone}} was run with 
#'   \code{eval=FALSE}, \code{get_design(x, 1)} will return the first
#'   normalization in the order saved by scone.
#'   
#' @export
#' 
setMethod(
  f = "get_design",
  signature = signature(x = "SconeExperiment", method = "numeric"),
  definition =  function(x, method) {
    norm_method <- rownames(x@scone_params)[method]

    return(get_design(x, norm_method))
  }
)
