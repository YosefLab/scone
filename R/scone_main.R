#' SCONE Main: Normalize Expression Data and Evaluate Normalization Performance
#'
#' This function is a high-level wrapper function for applying and evaluating
#' a variety of normalization schemes to a specified expression matrix.
#'
#' Each normalization consists of three main steps:
#' \itemize{
#'  \item{Impute:}{ Replace observations of zeroes with expected expression values. }
#'  \item{Scale:}{ Match sample-specific expression scales or quantiles. }
#'  \item{Adjust:}{ Adjust for sample-level batch factors / unwanted variation. }
#' }
#' Following completion of each step, the normalized expression matrix is scored
#' based on SCONE's data-driven evaluation criteria.
#' 
#' @param expr matrix. The expression data matrix (genes in rows, cells in columns).
#' @param imputation list or function. (A list of) function(s) to be used for
#'   imputation. By default only scone::impute_null is included.
#' @param impute_args arguments passed to all imputation functions.
#' @param rezero Restore zeroes following scaling step? Default FALSE.
#' @param scaling list or function. (A list of) function(s) to be used for
#'   scaling normalization.
#' @param k_ruv numeric. The maximum number of factors of unwanted variation
#'   (the function will adjust for 1 to k_ruv factors of unwanted variation). If
#'   0, RUV will not be performed.
#' @param k_qc numeric. The maximum number of quality metric PCs (the function
#'   will adjust for 1 to k_qc PCs). If 0, QC adjustment will not be performed.
#' @param ruv_negcon character. The genes to be used as negative controls for
#'   RUV. Ignored if k_ruv=0.
#' @param qc matrix. The QC metrics to be used for QC adjustment. Ignored if
#'   k_qc=0.
#' @param adjust_bio character. If 'no' it will not be included in the model; if
#'   'yes', both models with and without 'bio' will be run; if 'force', only
#'   models with 'bio' will be run.
#' @param adjust_batch character. If 'no' it will not be included in the model;
#'   if 'yes', both models with and without 'batch' will be run; if 'force',
#'   only models with 'batch' will be run.
#' @param bio factor. The biological condition to be included in the adjustment
#'   model (variation to be preserved). If adjust_bio="no", it will not be used
#'   for normalization, but only for evaluation.
#' @param batch factor. The known batch variable to be included in the
#'   adjustment model (variation to be removed). If adjust_batch="no", it will
#'   not be used for normalization, but only for evaluation.
#' @param evaluate logical. If FALSE the normalization methods will not be
#'   evaluated (faster).
#' @param eval_pcs numeric. The number of principal components to use for
#'   evaluation. Ignored if evaluation=FALSE.
#' @param eval_proj function. Projection function for evaluation  (see
#'   \code{\link{score_matrix}} for details). If NULL, PCA is used for
#'   projection.
#' @param eval_proj_args list. List of args passed to projection function as
#'   eval_proj_args.
#' @param eval_kclust numeric. The number of clusters (> 1) to be used for pam
#'   tightness evaluation. If an array of integers, largest average silhouette
#'   width (tightness) will be reported. If NULL, tightness will be returned NA.
#' @param eval_negcon character. The genes to be used as negative controls for
#'   evaluation. These genes should be expected not to change according to the
#'   biological phenomenon of interest. Ignored if evaluation=FALSE. Default is
#'   ruv_negcon argument. If NULL, correlations with negative controls will be
#'   returned NA.
#' @param eval_poscon character. The genes to be used as positive controls for
#'   evaluation. These genes should be expected to change according to the
#'   biological phenomenon of interest. Ignored if evaluation=FALSE. If NULL,
#'   correlations with positive controls will be returned NA.
#' @param params matrix or data.frame. If given, the algorithm will bypass
#'   creating the matrix of possible parameters, and will use the given matrix.
#'   There are basically no checks as to whether this matrix is in the right
#'   format, and is only intended to be used to feed the results of setting
#'   run=FALSE back into the algorithm (see example).
#' @param verbose logical. If TRUE some messagges are printed.
#' @param stratified_pam logical. If TRUE then maximum ASW for PAM_SIL is
#'   separately computed for each biological-cross-batch stratum (accepting
#'   NAs), and a weighted average is returned as PAM_SIL.
#' @param run logical. If FALSE the normalization and evaluation are not run,
#'   but the function returns a data.frame of parameters that will be run for
#'   inspection by the user.
#' @param return_norm character. If "no" the normalized values will not be
#'   returned. This will create a much smaller object and may be useful for
#'   large datasets and/or when many combinations are compared. If "in_memory"
#'   the normalized values will be returned as part of the output. If "hdf5"
#'   they will be written on file using the \code{rhdf5} package.
#' @param hdf5file character. If \code{return_norm="hdf5"}, the name of the file
#'   onto which to save the normalized matrices.
#' @param ... Arguments to be passed to methods.
#'
#' @return A list with the following elements: \itemize{ \item{normalized_data}{
#'   A list containing the normalized data matrix, log-scaled. NULL when
#'   \code{return_norm} is not "in_memory".} \item{metrics}{ A matrix containing
#'   raw evaluation metrics for each normalization method (see
#'   \code{\link{score_matrix}} for details). Rows are sorted in the same order
#'   as in the scores output matrix. NULL when evaluate = FALSE.} \item{scores}{
#'   A matrix containing scores for each normalization, including average score rank.
#'   Rows are sorted by decreasing mean score rank. NULL when evaluate = FALSE.}
#'   \item{params}{ A data.frame with each row corresponding to a set of
#'   normalization parameters.} }
#'   
#' @return If \code{run=FALSE} a \code{data.frame} with each row corresponding
#'   to a set of normalization parameters.
#'
#' @details Evaluation metrics are defined in \code{\link{score_matrix}}.
#'   Each metric is assigned a signature for conversion to scores:
#'   Positive-signature metrics increase with improving performance, including
#'   BIO_SIL, PAM_SIL, and EXP_WV_COR. Negative-signature metrics decrease with
#'   improving performance, including BATCH_SIL, EXP_QC_COR, EXP_UV_COR,
#'   RLE_MED, and RLE_IQR. Scores are computed so that higer-performing methods
#'   are assigned higher scores.
#'   
#' @details Note that if one wants to include the unnormalized data in the final
#'   comparison of normalized matrices, the identity function must be included 
#'   in the scaling list argument. Analogously, if one wants to include non-
#'   imputed data in the comparison, the scone::impute_null function must be 
#'   included.
#'
#' @details If \code{run=FALSE} the normalization and evaluation are not run,
#'   but the function returns a matrix of parameters that will be run for
#'   inspection by the user.
#'
#' @details If \code{return_norm="hdf5"}, the normalized matrices will be
#'   written to the \code{hdf5file} file. This must be a string specifying (a
#'   path to) a new file. If the file already exists, it will return error.
#'   
#' @importFrom RUVSeq RUVg
#' @importFrom matrixStats rowMedians
#' @import BiocParallel
#' @importFrom graphics abline arrows barplot hist par plot text
#' @importFrom stats approx as.formula binomial coefficients contr.sum cor dist
#'   dnbinom fitted.values glm mad median model.matrix na.omit p.adjust pnorm
#'   prcomp quantile quasibinomial sd
#' @importFrom utils capture.output
#' @importFrom rhdf5 h5createFile h5write.default h5write
#' @export
#' 

scone <- function(expr, imputation=list(none=impute_null), impute_args = NULL,
                  rezero = FALSE, scaling, k_ruv=5, k_qc=5,
                  ruv_negcon=NULL, qc=NULL, adjust_bio=c("no", "yes", "force"),
                  adjust_batch=c("no", "yes", "force"), bio=NULL, batch=NULL,
                  run=TRUE, evaluate=TRUE, eval_pcs=3, eval_proj = NULL,
                  eval_proj_args = NULL, eval_kclust=2:10, eval_negcon=ruv_negcon,
                  eval_poscon=NULL, params=NULL, verbose=FALSE,
                  stratified_pam = FALSE, return_norm = c("no", "in_memory", "hdf5"),
                  hdf5file, ...) {

  return_norm <- match.arg(return_norm)

  if(return_norm == "hdf5") {
    if(missing(hdf5file)) {
      stop("If `return_norm='hdf5'`, `hdf5file` must be specified.")
    } else {
      if(BiocParallel::bpnworkers(BiocParallel::registered()[[1]]) > 1) {
        stop("At the moment, `return_norm='hdf5'` does not support multicores.")
      }
      stopifnot(h5createFile(hdf5file))
      h5write(rownames(expr), hdf5file, "genes")
      h5write(colnames(expr), hdf5file, "samples")
    }
  }

  if(!is.matrix(expr)) {
    stop("'expr' must be a matrix.")
  } else if(is.null(rownames(expr))) {
    stop("'expr' must have row names.")
  }

  if(!is.function(imputation)) {
    if(is.list(imputation)) {
      if(!all(sapply(imputation, is.function))) {
        stop("'imputation' must be a function or a list of functions.")
      }
      if(is.null(names(imputation))) {
        names(imputation) <- paste("imputation", seq_along(imputation), sep="")
      }
    } else {
      stop("'imputation' must be a function or a list of functions.")
    }
  }

  if(is.function(imputation)) {
    l <- list(imputation)
    names(l) <- deparse(substitute(imputation))
    imputation <- l
  }

  if(!is.function(scaling)) {
    if(is.list(scaling)) {
      if(!all(sapply(scaling, is.function))) {
        stop("'scaling' must be a function or a list of functions.")
      }
      if(is.null(names(scaling))) {
        names(scaling) <- paste("scaling", seq_along(scaling), sep="")
      }
    } else {
      stop("'scaling' must be a function or a list of functions.")
    }
  }

  if(is.function(scaling)) {
    l <- list(scaling)
    names(l) <- deparse(substitute(scaling))
    scaling <- l
  }

  if(k_ruv < 0) stop("'k_ruv' must be non-negative.")
  if(k_qc < 0) stop("'k_qc' must be non-negative.")

  if(k_ruv > 0) {
    if(is.null(ruv_negcon)) {
      stop("If k_ruv>0, ruv_negcon must be specified.")
    } else if(!is.character(ruv_negcon)) {
      stop("'ruv_negcon' must be a character vector.")
    } else if(!all(ruv_negcon %in% rownames(expr))) {
      stop("'ruv_negcon' must be a subset of the genes in 'expr.'")
    }
  }

  if(k_qc > 0) {
    if(is.null(qc)) {
      stop("If k_qc>0, qc must be specified.")
    } else if(!is.matrix(qc)) {
      stop("'qc' must be a matrix.")
    } else if(nrow(qc) != ncol(expr)) {
      stop("'qc' must have one row per sample.")
    } else if(ncol(qc) < k_qc) {
      stop("k_qc must be less or equal than the number of columns of 'qc.'")
    }
    qc_pcs <- prcomp(qc, center=TRUE, scale=TRUE)$x
  } else {
    qc_pcs <- NULL
  }

  adjust_batch <- match.arg(adjust_batch)
  adjust_bio <- match.arg(adjust_bio)

  if(adjust_bio != "no") {
    if(is.null(bio)) {
      stop("if adjust_bio is 'yes' or 'force', 'bio' must be specified.")
    } else if(!is.factor(bio)) {
      stop("'bio' must be a factor.")
    } else if(length(bio) != ncol(expr)) {
      stop("'bio' must have length equal to the number of samples.")
    } else {
      bio <- factor(bio, levels = sort(levels(bio)))
    }
  }

  if(adjust_batch != "no") {
    if(is.null(batch)) {
      stop("if adjust_batch is 'yes' or 'force', 'batch' must be specified.")
    } else if(!is.factor(batch)) {
      stop("'batch' must be a factor.")
    } else if(length(batch) != ncol(expr)) {
      stop("'batch' must have length equal to the number of samples.")
    } else {
      if(!is.null(bio)) {
        batch <- factor(batch, levels = unique(batch[order(bio)]))
      } else {
        batch <- factor(batch, levels = sort(levels(batch)))
      }
    }
  }


  if(evaluate) {
    if(is.null(eval_negcon)) {
      if(verbose) message("eval_negcon is null, negative controls will not be used in evaluation (correlations with negative controls will be returned as NA)")
    } else if(!is.character(eval_negcon)) {
      stop("'eval_negcon' must be a character vector.")
    } else if(!all(eval_negcon %in% rownames(expr))) {
      stop("'eval_negcon' must be a subset of the genes in 'expr.'")
    }

    if(!is.character(eval_poscon) & !is.null(eval_poscon)) {
      stop("'eval_poscon' must be a character vector.")
    } else if(!all(eval_poscon %in% rownames(expr))) {
      stop("'eval_poscon' must be a subset of the genes in 'expr.'")
    }

    if(eval_pcs > ncol(expr)) {
      stop("'eval_pcs' must be less or equal than the number of samples.")
    }

    if(any(eval_kclust >= ncol(expr))) {
      stop("'eval_kclust' must be less than the number of samples.")
    }

    if(!is.null(eval_kclust) & stratified_pam) {
      if(is.null(bio) & is.null(batch)){
        stop("For stratified_pam, bio and/or batch must be specified")
      }
      biobatch = as.factor(paste(bio,batch,sep = "_"))
      if(max(eval_kclust) >= min(table(biobatch))) {
        stop("For stratified_pam, max 'eval_kclust' must be smaller than bio-cross-batch stratum size")
      }
    }
  }

  # Check Design: Confounded design or Nesting of Batch Condition and Biological Condition
  nested <- FALSE

  if(!is.null(bio) & !is.null(batch)) {
    tab <- table(bio, batch)
    if(all(colSums(tab>0)==1)){ # On
      if(nlevels(bio) == nlevels(batch)) {
        if(adjust_bio != "no" & adjust_batch != "no") {
          stop("Biological conditions and batches are confounded. They cannot both be included in the model, please set at least one of 'adjust_bio' and 'adjust_batch' to 'no.'")
        } else {
          warning("Biological conditions and batches are confounded. Removing batch effects may remove true biological signal and/or inferred differences may be inflated because of batch effects.")
        }
      } else {
        nested <- TRUE
        if(verbose) message("Detected nested design...")
      }
    }
  }

  # Step 0: compute the parameter matrix
  if(is.null(params)) {
    bi <- switch(adjust_bio,
                 no="no_bio",
                 yes=c("no_bio", "bio"),
                 force="bio"
    )

    ba <- switch(adjust_batch,
                 no="no_batch",
                 yes=c("no_batch", "batch"),
                 force="batch"
    )

    ruv_qc <- "no_uv"
    if(k_ruv > 0) {
      ruv_qc <- c(ruv_qc, paste("ruv_k=", 1:k_ruv, sep=""))
    }
    if(k_qc > 0) {
      ruv_qc <- c(ruv_qc, paste("qc_k=", 1:k_qc, sep=""))
    }

    params <- expand.grid(imputation_method=names(imputation),
                          scaling_method=names(scaling),
                          uv_factors=ruv_qc,
                          adjust_biology=bi,
                          adjust_batch=ba,
                          stringsAsFactors=FALSE
                          )
    rownames(params) <- apply(params, 1, paste, collapse=',')

    ## if adjust_bio = "yes", meaning both bio and no_bio, than remove bio when
    ## no batch or uv correction
    if(adjust_bio == "yes") {
      remove_params <- which(params$uv_factors=="no_uv" & params$adjust_batch=="no_batch"
                             & params$adjust_biology=="bio")
      params <- params[-remove_params,]
    }
  }

  if(!run) {
    return(params)
  }

  ## add a check to make sure that design matrix is full rank

  ## Step 1: imputation
  if(verbose) message("Imputation step...")
  im_params <- unique(params[,1])

  imputed <- lapply(seq_along(im_params), function(i) imputation[[im_params[i]]](expr,impute_args))
  names(imputed) <- im_params
  # output: a list of imputed matrices

  ## Step 2: scaling
  if(verbose) message("Scaling step...")
  sc_params <- unique(params[,1:2])

  scaled <- bplapply(seq_len(NROW(sc_params)), function(i) scaling[[sc_params[i,2]]](imputed[[sc_params[i,1]]]))
  if(rezero){
    if(verbose) message("Re-zero step...")
    toz = expr <= 0
    scaled <- lapply(scaled,function(x) x - x*toz)
  }
  names(scaled) <- apply(sc_params, 1, paste, collapse="_")
  # output: a list of normalized expression matrices

  failing <- bplapply(scaled, function(x) any(is.na(x)))
  failing <- simplify2array(failing)
  if(any(failing)) {
    idx <- which(failing)
    failed_names = names(scaled)[idx]
    warning(paste(failed_names, "returned at least one NA value. Consider removing it from the comparison.\n", 
                  "In the meantime, failed methods have been filtered from the output."))
    
    remove_params = paste(params$imputation_method, params$scaling_method, sep="_") %in% failed_names
    params <- params[!remove_params,]
    
    scaled = scaled[!failing]
  }

  if(verbose) message("Computing RUV factors...")
  ## compute RUV factors
  if(k_ruv > 0) {
    ruv_factors <- bplapply(scaled, function(x) RUVg(log1p(x), ruv_negcon, k_ruv, isLog=TRUE)$W)
  }

  if(evaluate) {

    if(verbose) message("Computing factors for evaluation...")

    ## generate factors: eval_pcs pcs per gene set
    uv_factors <- wv_factors <- NULL

    if(!is.null(eval_negcon)) {
      uv_factors <- svd(scale(t(log1p(expr[eval_negcon,])),center = TRUE,scale = TRUE),eval_pcs,0)$u
    }

    if(!is.null(eval_poscon)) {
      wv_factors <- svd(scale(t(log1p(expr[eval_poscon,])),center = TRUE,scale = TRUE),eval_pcs,0)$u
    }

  }

  if(verbose) message("Factor adjustment and evaluation...")

  outlist <- bplapply(seq_len(NROW(params)), function(i) {
    parsed <- parse_row(params[i,], bio, batch, ruv_factors, qc_pcs)
    design_mat <- make_design(parsed$bio, parsed$batch, parsed$W,
                              nested=(nested & !is.null(parsed$bio) & !is.null(parsed$batch)))
    sc_name <- paste(params[i,1:2], collapse="_")
    adjusted <- lm_adjust(log1p(scaled[[sc_name]]), design_mat, batch)
    if(evaluate) {
      #cat(sprintf("scoring matrix: I am in method: %s", sc_name)) #debug message
      score <- score_matrix(expr=adjusted, eval_pcs = eval_pcs, eval_proj = eval_proj, eval_proj_args = eval_proj_args,
                            eval_kclust = eval_kclust, bio = bio, batch = batch,
                            qc_factors = qc_pcs, uv_factors = uv_factors, wv_factors = wv_factors,
                            is_log = TRUE, stratified_pam = stratified_pam)
    } else {
      score <- NULL
    }

    if(verbose) message(paste0("Processed: ", rownames(params)[i]))
    if(return_norm == "in_memory") {
      return(list(score=score, adjusted=adjusted))
    } else {
      if(return_norm == "hdf5") {
        h5write(adjusted, hdf5file, rownames(params)[i])
      }
      return(list(score=score))
    }
  })

  if(return_norm == "in_memory") {
    adjusted <- lapply(outlist, function(x) x$adjusted)
    names(adjusted) <- apply(params, 1, paste, collapse=',')
  } else {
    adjusted <- NULL
  }

  if(evaluate) {
    evaluation <- lapply(outlist, function(x) x$score)
    names(evaluation) <- apply(params, 1, paste, collapse=',')
    evaluation <- simplify2array(evaluation)

    scores <- evaluation * c(1, -1, 1,  # "BIO_SIL", "BATCH_SIL", "PAM_SIL"
                             -1, -1, 1, # "EXP_QC_COR", "EXP_UV_COR", "EXP_WV_COR"
                             -1, -1) # "RLE_MED", "RLE_IQR"

    # Mean Score Rank per Method
    ranked_scores = apply(na.omit(scores),1,rank)
    if(is.null(dim(ranked_scores))){
      mean_score_rank <- mean(ranked_scores)
    }else{
      mean_score_rank <- rowMeans(ranked_scores)
    }
    
    scores <- cbind(t(scores), mean_score_rank)[order(mean_score_rank, decreasing = TRUE),]

    evaluation <- t(evaluation[,order(mean_score_rank, decreasing = TRUE), drop=FALSE])

    if(return_norm == "in_memory") {
      adjusted <- adjusted[order(mean_score_rank, decreasing = TRUE)]
    }

    params <- params[order(mean_score_rank, decreasing = TRUE),]
  } else {
    evaluation <- scores <- NULL
  }

  if(verbose) message("Done!")

  return(list(normalized_data=adjusted, metrics=evaluation, scores=scores, params=params))
}

#' Retrieve normalized matrix
#'
#' Given a string representing a normalization scheme and an HDF5 file created by
#' scone, it will return a matrix of normalized counts (in log scale).
#'
#' @details The string must be of of the row.names of the \code{params}
#'   component of the scone output.
#'
#' @param normalization character. A string identifying the normalization scheme
#'   to be retrieved.
#' @param hdf5file character. The name of the file in which the data have been
#'   stored.
#'
#' @return A matrix of normalized counts in log-scale.
#'
#' @importFrom rhdf5 h5ls h5read
#' @export
get_normalized <- function(normalization, hdf5file) {

  if(!normalization %in% h5ls(hdf5file)$name) {
    stop(paste0(normalization, " does not match any object in ", hdf5file))
  } else {
    retval <- h5read(hdf5file, normalization)
    rownames(retval) <- h5read(hdf5file, "genes")
    colnames(retval) <- h5read(hdf5file, "samples")
    return(retval)
  }
}
