#' scone evaluation: function to evaluate one normalization scheme
#'
#' This function evaluates an expression matrix using SCONE criteria, producing a number of scores based on
#' projections of the normalized data, correlations, and RLE metrics.
#'
#' @details The eval_proj function argument must have 2 inputs: 
#' \itemize{
#' \item{e}{ matrix. log-transformed expression (genes in rows, cells in columns).}
#' \item{eval_proj_args}{ list. additional function arguments, e.g. prior data weights.}
#' }
#' and it must output a matrix representation of the original data (cells in rows, factors in columns).
#' 
#' @param expr matrix. The data matrix (genes in rows, cells in columns).
#' @param eval_pcs numeric. The number of principal components to use for evaluation.
#' Ignored if !is.null(eval_proj).
#' @param eval_proj function. Projection function for evaluation (see details).
#' If NULL, PCA is used for projection
#' @param eval_proj_args list. List of arguments passed to projection function as eval_proj_args (see details).
#' @param eval_kclust numeric. The number of clusters (> 1) to be used for pam tightness (PAM_SIL) evaluation.
#' If an array of integers, largest average silhouette width (tightness) will be reported in PAM_SIL. If NULL, PAM_SIL will be returned NA.
#' @param bio factor. A known biological condition (variation to be preserved), NA is allowed.
#' If NULL, condition ASW, BIO_SIL, will be returned NA.
#' @param batch factor. A known batch variable (variation to be removed), NA is allowed.
#' If NULL, batch ASW, BATCH_SIL, will be returned NA.
#' @param qc_factors Factors of unwanted variation derived from quality metrics.
#' If NULL, qc correlations, EXP_QC_COR, will be returned NA.
#' @param uv_factors Factors of unwanted variation derived from negative control genes (evaluation set).
#' If NULL, uv correlations, EXP_UV_COR, will be returned NA.
#' @param wv_factors Factors of wanted variation derived from positive control genes (evaluation set).
#' If NULL, wv correlations, EXP_WV_COR, will be returned NA.
#' @param is_log logical. If TRUE the expr matrix is already logged and log transformation will not be carried out prior to projection.
#' @param conditional_pam logical. If TRUE then maximum ASW for PAM_SIL is separately computed for each biological condition (including NA),
#' and a weighted average silhouette width is returned.
#'
#' @importFrom class knn
#' @importFrom fpc pamk
#' @importFrom cluster silhouette
#' @importFrom matrixStats rowMedians colMedians colIQRs
#'
#' @export
#'
#' @return A list with the following elements:
#' \itemize{
#' \item{BIO_SIL}{ Average silhouette width by biological condition.}
#' \item{BATCH_SIL}{ Average silhouette width by batch condition.}
#' \item{PAM_SIL}{ Maximum average silhouette width from pam clustering (see conditional_pam argument).}
#' \item{EXP_QC_COR}{ Maximum squared spearman correlation between pcs and quality factors.}
#' \item{EXP_UV_COR}{ Maximum squared spearman correlation between pcs and passive uv factors.}
#' \item{EXP_WV_COR}{ Maximum squared spearman correlation between pcs and passive wv factors.}
#' \item{RLE_MED}{ The mean squared median Relative Log Expression (RLE).}
#' \item{RLE_IQR}{ The mean inter-quartile range (IQR) of the RLE.}
#' }
#'

score_matrix <- function(expr, eval_pcs = 3, 
                         eval_proj = NULL, eval_proj_args = NULL,
                         eval_kclust = NULL,
                         bio = NULL, batch = NULL,
                         qc_factors = NULL,uv_factors = NULL, wv_factors = NULL, 
                         is_log=FALSE, conditional_pam = FALSE){

  if(any(is.na(expr) | is.infinite(expr) | is.nan(expr))){
    stop("NA/Inf/NaN Expression Values.")
  }

  if(!is_log) {
    expr <- log1p(expr)
  }

  if(is.null(eval_proj)){
      proj = svd(scale(t(expr),center = TRUE,scale = TRUE),eval_pcs,0)$u
  }else{
      proj = eval_proj(expr,eval_proj_args = eval_proj_args)
      eval_pcs = ncol(proj)
  }

  ## ------ Bio and Batch Tightness -----

    if( !is.null(bio) ) {
      if(!all(is.na(bio))) {
        BIO_SIL = summary(cluster::silhouette(as.numeric(na.omit(bio)),dist(proj[!is.na(bio),])))$avg.width
      } else {
        BIO_SIL = NA
        warning("bio is all NA!")
      }
    } else {
      BIO_SIL = NA
    }

    if(!is.null(batch)) {
      if(!all(is.na(batch))) {
        BATCH_SIL <- summary(cluster::silhouette(as.numeric(na.omit(batch)),dist(proj[!is.na(batch),])))$avg.width
      } else{
        BATCH_SIL <- NA
        warning("batch is all NA!")
      }
    } else {
      BATCH_SIL <- NA
    }

  ## ------ PAM Tightness -----

  if ( !is.null(eval_kclust) ){
    
    # "Conditional" PAM
    if(conditional_pam){

      PAM_SIL = 0

      # Max Average Sil Width per Biological Condition
      for (cond in levels(bio)){
        is_cond = which(bio == cond)
        cond_w = length(is_cond)
        if(cond_w > max(eval_kclust)){
          pamk_object = pamk(proj[is_cond,],krange = eval_kclust)
          PAM_SIL = PAM_SIL + cond_w*pamk_object$pamobject$silinfo$avg.width
        }else{
          stop("Number of clusters 'k' must be smaller than bio class size")
        }
      }

      # Max Average Sil Width for Unclassified Condition
      is_na = is.na(bio)
      if (any(is_na)){
        cond_w = sum(is_na)
        if(cond_w > max(eval_kclust)){
          pamk_object = pamk(proj[is_na,],krange = eval_kclust)
          PAM_SIL = PAM_SIL + cond_w*pamk_object$pamobject$silinfo$avg.width
        }else{
          stop("Number of clusters 'k' must be smaller than unclassified bio set size")
        }
      }
      
      PAM_SIL = PAM_SIL/length(bio)
      
    # Traditional PAM
    }else{
      pamk_object = pamk(proj,krange = eval_kclust)
      PAM_SIL = pamk_object$pamobject$silinfo$avg.width
    }

  }else{
    PAM_SIL = NA
  }

  ## ------ Correlation with Factors -----

  # Max cor with quality factors.
  if(!is.null(qc_factors)){
    EXP_QC_COR = max(cor(proj,qc_factors,method = "spearman")^2)
  }else{
    EXP_QC_COR = NA
  }

  # Max cor with UV factors.
  if(!is.null(uv_factors)){
    EXP_UV_COR = max(cor(proj,uv_factors,method = "spearman")^2)
  }else{
    EXP_UV_COR = NA
  }

  # Max cor with WV factors.
  if(!is.null(wv_factors)){
    EXP_WV_COR = max(cor(proj,wv_factors,method = "spearman")^2)
  }else{
    EXP_WV_COR = NA
  }

  ## ----- RLE Measures
  rle <- expr - rowMedians(expr)
  RLE_MED <- mean(colMedians(rle)^2)
  RLE_IQR <- mean(colIQRs(rle))

  scores = c(BIO_SIL, BATCH_SIL, PAM_SIL, 
             EXP_QC_COR, EXP_UV_COR, EXP_WV_COR, 
             RLE_MED, RLE_IQR)
  names(scores) = c("BIO_SIL", "BATCH_SIL", "PAM_SIL", 
                    "EXP_QC_COR", "EXP_UV_COR", "EXP_WV_COR",
                    "RLE_MED", "RLE_IQR")
  return(scores)
}
