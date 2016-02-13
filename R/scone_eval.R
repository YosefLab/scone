#' scone evaluation: function to evaluate one normalization scheme
#' 
#' This function evaluates an expression matrix using SCONE criteria, producing a number of scores based on 
#' weighted (or unweighted) projections of the normalized data.
#' 
#' @details None.
#'  
#' @param expr matrix. The data matrix (genes in rows, cells in columns).
#' @param eval_pcs numeric. The number of principal components to use for evaluation.
#' Ignored if !is.null(proj).
#' @param proj matrix. A numeric data matrix to be used as projection (cells in rows, coordinates in columns).
#' If NULL, weighted PCA is used for projection
#' @param weights matrix. A numeric data matrix to be used for weighted PCA (genes in rows, cells in columns).
#' If NULL, regular PCA is used for projection
#' @param seed numeric. Random seed, passed to bwpca. 
#' Ignored if is.null(weights) or !is.null(proj).
#' @param em.maxiter numeric. Maximum EM iterations, passed to bwpca. 
#' Ignored if is.null(weights) or !is.null(proj).
#' @param eval_knn numeric. The number of nearest neighbors to use for evaluation.
#' If NULL, all KNN concordances will be returned NA.
#' @param eval_kclust numeric. The number of clusters (> 1) to be used for pam stability evaluation. 
#' If an array of integers, largest average silhoutte width will be reported. If NULL, stability will be returned NA.
#' @param bio factor. The biological condition (variation to be preserved), NA is allowed.
#' If NULL, condition KNN concordance will be returned NA.
#' @param batch factor. The known batch variable (variation to be removed), NA is allowed.
#' If NULL, batch KNN concordance will be returned NA.
#' @param qc_factors Factors of unwanted variation derived from quality metrics.
#' If NULL, qc correlations will be returned NA.
#' @param ruv_factors Factors of unwanted variation derived from negative control genes (adjustable set).
#' If NULL, ruv correlations will be returned NA.
#' @param uv_factors Factors of unwanted variation derived from negative control genes (un-adjustable set).
#' If NULL, uv correlations will be returned NA.
#' @param wv_factors Factors of wanted variation derived from positive control genes (un-adjustable set).
#' If NULL, wv correlations will be returned NA.
#' @param is_log logical. If TRUE the expr matrix is already logged and log transformation will not be carried out.
#' @param conditional_pam logical. If TRUE then maximum ASW is separately computed for each biological condition (including NA), 
#' and a weighted average is returned.
#' 
#' @importFrom scde bwpca
#' @importFrom class knn
#' @importFrom fpc pamk
#' 
#' @export
#' 
#' @return A list with the following elements:
#' \itemize{
#' \item{KNN_BIO}{ K-NN concordance rate by biological condition.}
#' \item{KNN_BATCH}{ K-NN concordance rate by batch condition.}
#' \item{PAM_SIL}{ Maximum average silhoutte width from pam clustering.}
#' \item{EXP_QC_COR}{ Maximum squared spearman correlation between pcs and quality factors.}
#' \item{EXP_RUV_COR}{ Maximum squared spearman correlation between pcs and active uv factors.}
#' \item{EXP_UV_COR}{ Maximum squared spearman correlation between pcs and passive uv factors.}
#' \item{EXP_WV_COR}{ Maximum squared spearman correlation between pcs and passive wv factors.}
#' }
#'

score_matrix <- function(expr, eval_pcs = 3, proj = NULL, 
                        weights = NULL, seed = 1, em.maxiter = 100,
                        eval_knn = NULL, eval_kclust = NULL,
                        bio = NULL, batch = NULL,
                        qc_factors = NULL, 
                        ruv_factors = NULL, uv_factors = NULL, 
                        wv_factors = NULL, is_log=FALSE, conditional_pam = FALSE ){
  
  if(any(is.na(expr) | is.infinite(expr) | is.nan(expr))){
    stop("NA/Inf/NaN Expression Values.")
  }
  
  if(!is_log) {
    expr <- log1p(expr)
  }
  
  if(is.null(proj)){
    if(is.null(weights)){
      proj = prcomp(t(expr),center = TRUE,scale. = TRUE)$x[,1:eval_pcs]
    }else{
      proj = bwpca(mat = t(expr),matw = t(weights),npcs = eval_pcs, seed = seed, em.maxiter = em.maxiter)$scores
    }
  }else{
    eval_pcs = dim(proj)[2]
  }
  
  ## ------ K-nearest neighbors (including self!) -----
  
  if( !is.null(eval_knn)  ){
    
    if( !is.null(bio) | !any(!is.na(bio)) ){
      KNN_BIO = mean(attributes(knn(train = proj[!is.na(bio),],test = proj[!is.na(bio),],cl = bio[!is.na(bio)], k = eval_knn,prob = TRUE))$prob)
    }else{
      KNN_BIO = NA
    }
    
  
    # K-NN Batch
    if( !is.null(batch) | !any(!is.na(batch)) ){
      KNN_BATCH = mean(attributes(knn(train = proj[!is.na(batch),],test = proj[!is.na(batch),],cl = batch[!is.na(batch)], k = eval_knn,prob = TRUE))$prob)
    }else{
      KNN_BATCH = NA
    }
    
  }else{
    
    KNN_BIO = NA
    KNN_BATCH = NA
    
  }
  
  ## ------ PAM Stability -----
  
  if ( !is.null(eval_kclust) ){
    
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
    }else{
      pamk_object = pamk(proj,krange = eval_kclust)
      PAM_SIL = pamk_object$pamobject$silinfo$avg.width
    }
  }else{
    PAM_SIL = NA
  }
  
  ## ------ Hidden Factors -----
  
  # Max cor with quality factors.
  if(!is.null(qc_factors)){
    EXP_QC_COR = max(cor(proj,qc_factors,method = "spearman")^2)  
  }else{
    EXP_QC_COR = NA
  }
  
  # Max cor with RUV factors.
  if(!is.null(ruv_factors)){
    EXP_RUV_COR = max(cor(proj,ruv_factors,method = "spearman")^2)  
  }else{
    EXP_RUV_COR = NA
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
    
  scores = c(KNN_BIO, KNN_BATCH, PAM_SIL, EXP_QC_COR, EXP_RUV_COR, EXP_UV_COR, EXP_WV_COR)
  names(scores) = c("KNN_BIO", "KNN_BATCH", "PAM_SIL", "EXP_QC_COR", "EXP_RUV_COR", "EXP_UV_COR", "EXP_WV_COR")
  return(scores)
  
}
