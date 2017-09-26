#' SCONE Evaluation: Evaluate an Expression Matrix
#'
#' This function evaluates a (normalized) expression matrix using SCONE
#' criteria, producing 8 metrics based on i) Clustering, ii) Correlations and
#' iii) Relative Expression.
#'
#' @details Users may specify their own eval_proj function that will be used to
#'   compute Clustering and Correlation metrics. This eval_proj() function must
#'   have 2 input arguments: \itemize{ \item{e}{ matrix. log-transformed (+
#'   pseudocount) expression data (genes in rows, cells in columns).}
#'   \item{eval_proj_args}{ list. additional function arguments, e.g. prior
#'   data weights.}} and it must output a matrix representation of the original
#'   data (cells in rows, factors in columns). The value of eval_proj_args is
#'   passed to the user-defined function from the eval_proj_args argument of
#'   the main score_matrix() function call.
#'
#' @param expr matrix. The expression data matrix (genes in rows, cells in
#'   columns).
#' @param eval_pcs numeric. The number of principal components to use for
#'   evaluation (Default 3). Ignored if !is.null(eval_proj).
#' @param eval_proj function. Projection function for evaluation (see Details).
#'   If NULL, PCA is used for projection
#' @param eval_proj_args list. List of arguments passed to projection function
#'   as eval_proj_args (see Details).
#' @param eval_kclust numeric. The number of clusters (> 1) to be used for pam
#'   tightness (PAM_SIL) evaluation. If an array of integers, largest average
#'   silhouette width (tightness) will be reported in PAM_SIL. If NULL, PAM_SIL
#'   will be returned NA.
#' @param bio factor. A known biological condition (variation to be preserved),
#'   NA is allowed. If NULL, condition ASW, BIO_SIL, will be returned NA.
#' @param batch factor. A known batch variable (variation to be removed), NA is
#'   allowed. If NULL, batch ASW, BATCH_SIL, will be returned NA.
#' @param qc_factors Factors of unwanted variation derived from quality
#'   metrics. If NULL, qc correlations, EXP_QC_COR, will be returned NA.
#' @param uv_factors Factors of unwanted variation derived from negative
#'   control genes (evaluation set). If NULL, uv correlations, EXP_UV_COR,
#'   will be returned NA.
#' @param wv_factors Factors of wanted variation derived from positive
#'   control genes (evaluation set). If NULL, wv correlations, EXP_WV_COR,
#'   will be returned NA.
#' @param is_log logical. If TRUE the expr matrix is already logged and log
#'   transformation will not be carried out prior to projection. Default FALSE.
#' @param stratified_pam logical. If TRUE then maximum ASW is separately
#'   computed for each biological-cross-batch stratum (accepts NAs), and a
#'   weighted average silhouette width is returned as PAM_SIL. Default FALSE.
#' @param stratified_cor logical. If TRUE then cor metrics are separately
#'   computed for each biological-cross-batch stratum (accepts NAs), and
#'   weighted averages are returned for EXP_QC_COR, EXP_UV_COR, & EXP_WV_COR.
#'   Default FALSE.
#' @param stratified_rle logical. If TRUE then rle metrics are separately
#'   computed for each biological-cross-batch stratum (accepts NAs), and
#'   weighted averages are returned for RLE_MED & RLE_IQR. Default FALSE.
#'
#' @importFrom class knn
#' @importFrom fpc pamk
#' @importFrom cluster silhouette
#' @importFrom rARPACK svds
#' @importFrom matrixStats rowMedians colMedians colIQRs
#'
#' @export
#'
#' @return A list with the following metrics: \itemize{ \item{BIO_SIL}{ Average
#'   silhouette width by biological condition.} \item{BATCH_SIL}{ Average
#'   silhouette width by batch condition.} \item{PAM_SIL}{ Maximum average
#'   silhouette width from PAM clustering (see stratified_pam argument).}
#'   \item{EXP_QC_COR}{ Coefficient of determination between expression
#'   pcs and quality factors (see stratified_cor argument).} \item{EXP_UV_COR}{
#'   Coefficient of determination between expression pcs and negative
#'   control gene factors (see stratified_cor argument).} \item{EXP_WV_COR}{
#'   Coefficient of determination between expression pcs and positive
#'   control gene factors (see stratified_cor argument).} \item{RLE_MED}{ The
#'   mean squared median Relative Log Expression (RLE) (see stratified_rle
#'   argument).} \item{RLE_IQR}{ The variance of the inter-quartile range (IQR)
#'   of the RLE (see stratified_rle argument).} }
#'
#' @examples
#'
#' set.seed(141)
#' bio = as.factor(rep(c(1,2),each = 2))
#' batch = as.factor(rep(c(1,2),2))
#' log_expr = matrix(rnorm(20),ncol = 4)
#'
#' scone_metrics = score_matrix(log_expr,
#'    bio = bio, batch = batch,
#'    eval_kclust = 2, is_log = TRUE)
#'

score_matrix <- function(expr,
                         eval_pcs = 3,
                         eval_proj = NULL,
                         eval_proj_args = NULL,
                         eval_kclust = NULL,
                         bio = NULL,
                         batch = NULL,
                         qc_factors = NULL,
                         uv_factors = NULL,
                         wv_factors = NULL,
                         is_log = FALSE,
                         stratified_pam = FALSE,
                         stratified_cor = FALSE,
                         stratified_rle = FALSE) {
  if (any(is.na(expr) | is.infinite(expr) | is.nan(expr))) {
    stop("NA/Inf/NaN Expression Values.")
  }
  
  if (!is_log) {
    expr <- log1p(expr)
  }
  
  # The svd we do below on expr throws an exception if
  # expr created by one of the normalizations has a
  # constant feature (=gene, i.e. row)
  constantFeatures = apply(expr, 1, function(x)
    max(x) - min(x)) < 1e-3
  if (any(constantFeatures)) {
    warning(sprintf(
      paste0(
        "scone_eval: expression matrix ",
        "contained %d constant features (rows) ",
        "---> excluding them"
      ),
      sum(constantFeatures)
    ))
    expr = expr[!constantFeatures,]
  }
  
  if (is.null(eval_proj)) {
    proj = tryCatch({
      svds(
        scale(t(expr), center = TRUE, scale = TRUE),
        k = eval_pcs,
        nu = eval_pcs,
        nv = 0
      )$u
    },
    error = function(e) {
      stop("scone_eval: svd failed")
    })
    
    
  } else {
    proj = eval_proj(expr, eval_proj_args = eval_proj_args)
    eval_pcs = ncol(proj)
  }
  
  ## ------ Bio and Batch Tightness -----
  dd <- as.matrix(dist(proj))
  
  # Biological Condition
  
  if (!is.null(bio)) {
    if (!all(is.na(bio))) {
      if (length(unique(bio)) > 1) {
        BIO_SIL = summary(cluster::silhouette(as.numeric(na.omit(bio)),
                                              dd[!is.na(bio),
                                                 !is.na(bio)]))$avg.width
      } else {
        BIO_SIL = NA
        warning(
          paste0(
            "after exclusion of samples, ",
            "only one bio remains, BATCH_BIO",
            " is undefined"
          )
        )
      }
    } else {
      BIO_SIL = NA
      warning("bio is all NA!")
    }
  } else {
    BIO_SIL = NA
  }
  
  # Batch Condition
  
  if (!is.null(batch)) {
    if (!all(is.na(batch))) {
      if (length(unique(batch)) > 1) {
        BATCH_SIL <- summary(cluster::silhouette(as.numeric(na.omit(batch)),
                                                 dd[!is.na(batch),
                                                    !is.na(batch)]))$avg.width
      } else {
        BATCH_SIL <- NA
        warning(
          paste0(
            "after exclusion of samples,",
            " only one batch remains, ",
            "BATCH_SIL is undefined"
          )
        )
      }
    } else{
      BATCH_SIL <- NA
      warning("batch is all NA!")
    }
  } else {
    BATCH_SIL <- NA
  }
  
  ## ------ PAM Tightness -----
  
  if (!is.null(eval_kclust)) {
    # "Stratified" PAM
    
    if (stratified_pam) {
      biobatch = as.factor(paste(bio, batch, sep = "_"))
      PAM_SIL = 0
      
      # Max Average Sil Width per bio-cross-batch Condition
      for (cond in levels(biobatch)) {
        is_cond = which(biobatch == cond)
        cond_w = length(is_cond)
        if (cond_w > max(eval_kclust)) {
          pamk_object = pamk(proj[is_cond, ], krange = eval_kclust)
          
          # Despite krange excluding nc = 1,
          # if asw is negative, nc = 1 will be selected
          if (is.null(pamk_object$pamobject$silinfo$avg.width)) {
            if (!1 %in% eval_kclust) {
              pamk_object$pamobject$silinfo$avg.width = 
                max(pamk_object$crit[1:max(eval_kclust) %in% eval_kclust])
            } else{
              stop(paste0(
                "nc = 1 was selected by Duda-Hart,",
                " exclude 1 from eval_kclust."
              ))
            }
          }
          
          PAM_SIL = PAM_SIL + cond_w * pamk_object$pamobject$silinfo$avg.width
        } else{
          stop(paste(
            paste0(
              "Number of clusters 'k' must be ",
              "smaller than bio-cross-batch stratum size:"
            ),
            paste(
              levels(biobatch),
              table(biobatch),
              sep = " = ",
              collapse = ", "
            )
          ))
        }
      }
      PAM_SIL = PAM_SIL / length(biobatch)
      
      # Traditional PAM
      
    } else{
      pamk_object = pamk(proj, krange = eval_kclust)
      PAM_SIL = pamk_object$pamobject$silinfo$avg.width
    }
    
  } else{
    PAM_SIL = NA
  }
  
  ## ------ Correlation with Factors -----
  if (stratified_cor) {
    biobatch = as.factor(paste(bio, batch, sep = "_"))
    EXP_QC_COR <- EXP_UV_COR <- EXP_WV_COR <- 0
    
    for (cond in levels(biobatch)) {
      is_cond = which(biobatch == cond)
      cond_w = length(is_cond)
      
      # Max cor with quality factors.
      if (!is.null(qc_factors)) {
        EXP_QC_COR <-
          EXP_QC_COR + cond_w * (1 - sum(unlist(apply(proj[is_cond,], 2, function(y) {
            lm(y ~ qc_factors[is_cond,])$residual
          })) ^ 2) / sum(scale(proj[is_cond,], scale = FALSE) ^ 2))
      }
      
      # Max cor with UV factors.
      if (!is.null(uv_factors)) {
        EXP_UV_COR  <-
          EXP_UV_COR + cond_w * (1 - sum(unlist(apply(proj[is_cond,], 2, function(y) {
            lm(y ~ uv_factors[is_cond,])$residual
          })) ^ 2) / sum(scale(proj[is_cond,], scale = FALSE) ^ 2))
      }
      
      # Max cor with WV factors.
      if (!is.null(wv_factors)) {
        EXP_WV_COR <-
          EXP_WV_COR + cond_w * (1 - sum(unlist(apply(proj[is_cond,], 2, function(y) {
            lm(y ~ wv_factors[is_cond,])$residual
          })) ^ 2) / sum(scale(proj[is_cond,], scale = FALSE) ^ 2))
      }
      
    }
    
    if (!is.null(qc_factors)) {
      EXP_QC_COR <- EXP_QC_COR / length(biobatch)
    } else{
      EXP_QC_COR <- NA
    }
    if (!is.null(uv_factors)) {
      EXP_UV_COR <- EXP_UV_COR / length(biobatch)
    } else{
      EXP_UV_COR <- NA
    }
    if (!is.null(wv_factors)) {
      EXP_WV_COR <- EXP_WV_COR / length(biobatch)
    } else{
      EXP_WV_COR <- NA
    }
    
  } else{
    # Max cor with quality factors.
    if (!is.null(qc_factors)) {
      EXP_QC_COR <- 1 - sum(unlist(apply(proj, 2, function(y) {
        lm(y ~ qc_factors)$residual
        })) ^ 2) / sum(scale(proj, scale = FALSE) ^ 2)
    } else{
      EXP_QC_COR = NA
    }
    
    # Max cor with UV factors.
    if (!is.null(uv_factors)) {
      EXP_UV_COR <- 1 - sum(unlist(apply(proj, 2, function(y) {
        lm(y ~ uv_factors)$residual
      })) ^ 2) / sum(scale(proj, scale = FALSE) ^ 2)
    } else{
      EXP_UV_COR = NA
    }
    
    # Max cor with WV factors.
    if (!is.null(wv_factors)) {
      EXP_WV_COR <- 1 - sum(unlist(apply(proj, 2, function(y) {
        lm(y ~ wv_factors)$residual
      })) ^ 2) / sum(scale(proj, scale = FALSE) ^ 2)
    } else{
      EXP_WV_COR = NA
    }
    
  }
  
  ## ----- RLE Measures
  if (stratified_rle) {
    biobatch = as.factor(paste(bio, batch, sep = "_"))
    RLE_MED <- RLE_IQR <- 0
    
    for (cond in levels(biobatch)) {
      is_cond = which(biobatch == cond)
      cond_w = length(is_cond)
      rle <- expr[, is_cond] - rowMedians(expr[, is_cond])
      
      # Non-Zero Median RLE
      RLE_MED <-
        RLE_MED + cond_w * mean(colMedians(rle) ^ 2) # Var of the med
      
      # Variable IQR RLE
      RLE_IQR <-
        RLE_IQR + cond_w * var(colIQRs(rle)) # Variance of the IQR
      
    }
    
    RLE_MED <- RLE_MED / length(biobatch)
    RLE_IQR <- RLE_IQR / length(biobatch)
    
  } else{
    rle <- expr - rowMedians(expr)
    
    # Non-Zero Median RLE
    RLE_MED <- mean(colMedians(rle) ^ 2) # Variance of the median
    
    # Variable IQR RLE
    RLE_IQR <- var(colIQRs(rle)) # Variance of the IQR
  }
  
  scores = c(BIO_SIL,
             BATCH_SIL,
             PAM_SIL,
             EXP_QC_COR,
             EXP_UV_COR,
             EXP_WV_COR,
             RLE_MED,
             RLE_IQR)
  names(scores) = c(
    "BIO_SIL",
    "BATCH_SIL",
    "PAM_SIL",
    "EXP_QC_COR",
    "EXP_UV_COR",
    "EXP_WV_COR",
    "RLE_MED",
    "RLE_IQR"
  )
  return(scores)
}
