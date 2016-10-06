#' @describeIn get_normalized If \code{method} is a character, it will return
#'   the normalized matrix corresponding to the normalization scheme specified
#'   by the character string. The string must be one of the \code{row.names} of
#'   the slot \code{scone_params}.
#'
#' @importFrom rhdf5 h5ls h5read
#' @export
#'
#' @param log logical. Should the data be returned in log-scale
#'
setMethod(
  f = "get_normalized",
  signature = signature(x = "SconeExperiment", method = "character"),
  definition =  function(x, method, log=FALSE) {

    if(!(method %in% rownames(x@scone_params))) {
      stop("`method` must be one of the row.names of the slot scone_params.")
    }

    if(x@scone_run == "in_memory") {
      retval <- assay(x, method)
    }

    if(x@scone_run == "hdf5") {
      retval <- h5read(hdf5file, method)
      rownames(retval) <- h5read(hdf5file, "genes")
      colnames(retval) <- h5read(hdf5file, "samples")
    }

    if(x@scone_run == "no") {
      params <- unlist(x@scone_params[method,])

      imputed <- x@imputation_fn[params[1]](assay(x))
      scaled <- x@scaling_fn[params[2]](imputed)

      ruv_factors <- qc_factors <- NULL
      if(params[3] != "no_uv") {

        k <- as.numeric(strsplit(params[3], "=")[[1]][2])

        if(grepl("ruv", params[3])) {
          ruv_factors <- RUVg(log1p(scaled), get_negconruv(x), k, isLog=TRUE)$W
        }

        if(grepl("qc", params[3])) {
          qc_factors <- prcomp(get_qc(x), center=TRUE, scale=TRUE)$x[,seq_len(k)]
        }

      }

      parsed <- parse_row(params, get_bio(x), get_batch(x), ruv_factors, qc_factors)
      design_mat <- make_design(parsed$bio, parsed$batch, parsed$W,
                                nested=(x@nested & !is.null(parsed$bio) & !is.null(parsed$batch)))
      adjusted <- lm_adjust(log1p(scaled), design_mat, get_batch(x))
      retval <- exp1m(adjusted)
    }

    if(log) {
      retval <- log1p(retval)
    }

    return(retval)
  }
)


#' @describeIn get_normalized If \code{method} is a numeric, it will return the normalized matrix according to the scone ranking.
#'
#' @export
#'
setMethod(
  f = "get_normalized",
  signature = signature(x = "SconeExperiment", method = "numeric"),
  definition =  function(x, method, log=FALSE) {
    norm_method <- rownames(x@scone_params)[method]

    return(get_normalized(x, norm_method, log))
  }
)

