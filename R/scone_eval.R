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
#' @param eval_kclust numeric. The number of clusters (> 1) to be used for pam tightness and stability evaluation.
#' If an array of integers, largest average silhoutte width (tightness) / maximum co-clustering compactness (stability) will be reported. If NULL, tightness and stability will be returned NA.
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
#' @param ref_expr matrix. A reference (log-) expression data matrix for calculating preserved variance (genes in rows, cells in columns).
#' If NULL, preserved variance is returned NA.
#'  
#' @importFrom scde bwpca
#' @importFrom class knn
#' @importFrom fpc pamk
#' @importFrom clusterCells subsampleClustering
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
#' \item{PAM_COMPACT}{ Compactness measure of sub-sampled (pam) co-clustering matrix's "block-diagonal-ness". Approximate isoperimetric quotient of non-clustering region.}
#' \item{VAR_PRES}{ Variance preserved measure.}
#' }
#'

score_matrix <- function(expr, eval_pcs = 3, proj = NULL,
                        weights = NULL, seed = 1, em.maxiter = 100,
                        eval_knn = NULL, eval_kclust = NULL,
                        bio = NULL, batch = NULL,
                        qc_factors = NULL,
                        ruv_factors = NULL, uv_factors = NULL,
                        wv_factors = NULL, is_log=FALSE, 
                        conditional_pam = FALSE , cv_genes = NULL, ref_expr = NULL){

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

    if( !is.null(bio) ) {
      if(!all(is.na(bio))) {
        KNN_BIO = mean(attributes(knn(train = proj[!is.na(bio),],test = proj[!is.na(bio),],cl = bio[!is.na(bio)], k = eval_knn,prob = TRUE))$prob)
      } else {
        KNN_BIO = NA
        warning("bio is all NA!")
      }
    } else {
      KNN_BIO = NA
    }

    # K-NN Batch
    if(!is.null(batch)) {
      if(!all(is.na(batch))) {
        KNN_BATCH <- mean(attributes(knn(train = proj[!is.na(batch),],test = proj[!is.na(batch),],cl = batch[!is.na(batch)], k = eval_knn,prob = TRUE))$prob)
      } else{
        KNN_BATCH <- NA
        warning("batch is all NA!")
      }
    } else {
      KNN_BATCH <- NA
    }

  } else {

    KNN_BIO = NA
    KNN_BATCH = NA

  }

  ## ------ PAM Tightness and Stability -----

  if ( !is.null(eval_kclust) ){

    # Tightness

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

    # Stability

    PAM_COMPACT = 0
    for(k in eval_kclust){

      submat = subsampleClustering(proj, k=k) # Co-clustering frequency matrix
      subhc = hclust(dist(submat)) # Re-order rows and columns using hclust
      submat = submat[subhc$order,subhc$order]

      submat_shifted = 2*submat - 1
      field = (submat_shifted[2:(nrow(submat_shifted)-1),3:nrow(submat_shifted)] + submat_shifted[2:(nrow(submat_shifted)-1),1:(nrow(submat_shifted)-2)])/2
      spin = submat_shifted[2:(nrow(submat_shifted)-1),2:(nrow(submat_shifted)-1)]
      perim_len = mean(1/2-field*spin/2) # Approximate fraction of elements on the perimeter
      isoper = (mean(1-submat)/(perim_len^2))/length(submat) # Approximate isoperimetric quotient of non-clustering (0) region.
      PAM_COMPACT = max(PAM_COMPACT,isoper) # Maximum compactness across all choices of resampling scheme.

    }

  }else{
    PAM_SIL = NA
    PAM_COMPACT = NA
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
  
  # VAR_PRES
  if(!is.null(ref_expr)){
    z1 = scale(t(ref_expr),scale = FALSE)
    z1 = t(t(z1)/sqrt(colSums(z1^2)))
    z2 = scale(t(expr),scale = FALSE)
    z2 = t(t(z2)/sqrt(colSums(z2^2)))
    VAR_PRES = mean(diag(t(z1) %*% (z2)))
  }else{
    VAR_PRES = NA
  }

  scores = c(KNN_BIO, KNN_BATCH, PAM_SIL, EXP_QC_COR, EXP_RUV_COR, EXP_UV_COR, EXP_WV_COR , PAM_COMPACT, VAR_PRES)
  names(scores) = c("KNN_BIO", "KNN_BATCH", "PAM_SIL", "EXP_QC_COR", "EXP_RUV_COR", "EXP_UV_COR", "EXP_WV_COR", "PAM_COMPACT","VAR_PRES")
  return(scores)

}
