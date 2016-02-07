#' Fit Logistic Regression Model of FNR against set of positive control (ubiquitously expressed) genes
#' 
#' @details logit(Probability of False Negative) ~ a + b*(mean log10p1 expression) .
#'  
#' @param expr matrix. The data matrix (genes in rows, cells in columns).
#' @param pos_controls. A boolean vector indexing positive control genes that will be used to compute false-negative rate characteristics.
#' @param fn_tresh. Inclusive threshold for negative detection. Default 0.01.
#' 
#' @return A list of logistic regression coefficients corresponding to glm fits in each sample. If a fit did not converge, the result reported is NA.
#'

simple_FNR_params = function(expr, pos_controls, fn_tresh = 0.01){
  
  # Mean log10p1 expression
  mu_obs = rowMeans(log10(expr[pos_controls,])+1)
  
  # Drop-outs
  drop_outs = 0 + (expr[pos_controls,] <= fn_tresh)
  
  # Logistic Regression Model of FNR
  ref.glms = list()
  for (si in 1:dim(drop_outs)[2]){
    fit = suppressWarnings(glm(cbind(drop_outs[,si],1 - drop_outs[,si]) ~ mu_obs,family=binomial(logit)))
    if(fit$converged){
      ref.glms[[si]] = fit$coefficients
    } else {
      ref.glms[[si]] = NA
    }
  }
  return(ref.glms)
}

#' metric-based sample filtering: function to filter single-cell RNA-Seq libraries.
#' 
#' This function returns a sample-filtering report for each cell in the input expression matrix, 
#' describing which filtering criteria are satisfied.
#' 
#' @details For each primary criterion (metric), a sample is evaluated based on 4 sub-criteria: 
#' 1) Hard (encoded) threshold 
#' 2) Adaptive thresholding via sd's from the mean
#' 3) Adaptive thresholding via mad's from the median
#' 4) Adaptive thresholding via sd's from the mean (after mixture modeling)
#' A sample must pass all sub-criteria to pass the primary criterion.
#'  
#' @param expr matrix. The data matrix (genes in rows, cells in columns).
#' @param nreads. A numeric vector representing number of reads in each library.
#' If NULL, filtered_nreads will be returned NA.
#' @param ralign. A numeric vector representing the proportion of reads aligned to the reference genome in each library.
#' If NULL, filtered_ralign will be returned NA.
#' @param gene_filter. A boolean vector indexing genes that will be used to compute library transcriptome breadth.
#' If NULL, filtered_breadth will be returned NA.
#' @param pos_controls. A boolean vector indexing positive control genes that will be used to compute false-negative rate characteristics.
#' If NULL, filtered_fnr will be returned NA.
#' @param AUC_range. An array of two values, representing range over which FNR AUC will be computed (log10(expr_units + 1)). Default c(0,6)
#' @param zcut. A numeric value determining threshold Z-score for sd, mad, and mixture sub-criteria. Default 1.
#' If NULL, only hard threshold sub-criteria will be applied.
#' @param mixture. A boolean value determining whether mixture modeling sub-criterion will be applied per primary criterion (metric).
#' If true, a dip test will be applied to each metric. If a metric is multimodal, it is fit to a two-component nomal mixture model. 
#' Samples deviating zcut sd's from optimal mean (in the inferior direction), have failed this sub-criterion.
#' @param dip_thresh. A numeric value determining dip test p-value threshold. Default 0.05.
#' @param hard_nreads numeric. Hard (lower bound on) nreads threshold. Default 25000.
#' @param hard_ralign numeric. Hard (lower bound on) ralign threshold. Default 15.
#' @param hard_breadth numeric. Hard (lower bound on) breadth threshold. Default 0.2.
#' @param hard_fnr numeric. Hard (upper bound on) fnr threshold. Default 3.
#' @param suff_nreads numeric. If not null, serves as an upper bound on nreads threshold.
#' @param suff_ralign numeric. If not null, serves as an upper bound on ralign threshold.  Default 65.
#' @param suff_breadth numeric. If not null, serves as an upper bound on breadth threshold.  Default 0.8.
#' @param suff_fnr numeric. If not null, serves as an lower bound on fnr threshold.
#' 
#' @return A list with the following elements:
#' \itemize{
#' \item{filtered_nreads}{Boolean. Sample has too few reads.}
#' \item{filtered_ralign}{Boolean. Sample has too few reads aligned.}
#' \item{filtered_breadth}{Boolean. Samples has too few genes detected (low breadth).}
#' \item{filtered_fnr}{Boolean. Sample has a high FNR AUC.}
#' }
#'

metric_sample_filter = function(expr, nreads = NULL, ralign = NULL,
                                gene_filter = NULL, pos_controls = NULL,
                                AUC_range = c(0,6), zcut = 1,
                                mixture = TRUE, dip_thresh = 0.05, 
                                hard_nreads = 25000, hard_ralign = 15, hard_breadth = 0.2, hard_fnr = 3,
                                suff_nreads = NULL, suff_ralign = 65, suff_breadth = 0.8, suff_fnr = NULL){
    
  ### ----- Primary Criterion 1) Number of Reads. -----
  
  if(!is.null(nreads)){ 
        
    logr = log10(nreads + 1)
  
    LOGR_CUTOFF = log10(hard_nreads + 1) # Hard Sub-Criterion
    
    if (!is.null(zcut)){
      LOGR_CUTOFF = max(mean(logr) - zcut*sd(logr), LOGR_CUTOFF) # SD Sub-Criterion
      LOGR_CUTOFF = max(median(logr) - zcut*mad(logr), LOGR_CUTOFF ) # MAD Sub-Criterion
    
      if(mixture){  # Mixture SD Sub-Criterion

        is.multimodal = dip.test(logr)$p.value < dip_thresh
        if(is.multimodal){
          mixmdl = normalmixEM(logr,k=2)
          component = which(mixmdl$mu %in% max(mixmdl$mu))
          LOGR_CUTOFF = max(mixmdl$mu[component] - zcut*mixmdl$sigma[component], LOGR_CUTOFF)
        }
      }

      if(!is.null(suff_nreads)){
        LOGR_CUTOFF = min(LOGR_CUTOFF,log10(suff_nreads + 1))
      }
    }
  }
  
  ## ----- Primary Criterion 2) Ratio of reads aligned. -----
  
  if(!is.null(ralign)){
    
    RALIGN_CUTOFF = hard_ralign
    
    if (!is.null(zcut)){
      RALIGN_CUTOFF = max(mean(ralign) - zcut*sd(ralign), RALIGN_CUTOFF) # SD Sub-Criterion
      RALIGN_CUTOFF = max(median(ralign) - zcut*mad(ralign), RALIGN_CUTOFF) # MAD Sub-Criterion
        
      if(mixture){ # Mixture SD Sub-Criterion
        
        is.multimodal = dip.test(ralign)$p.value < dip_thresh
        
        if(is.multimodal){
          mixmdl = normalmixEM(ralign,k=2)
          component = which(mixmdl$mu %in% max(mixmdl$mu))
          RALIGN_CUTOFF = max(mixmdl$mu[component] - zcut*mixmdl$sigma[component], RALIGN_CUTOFF)
        }
      }
      
      if(!is.null(suff_ralign)){
        RALIGN_CUTOFF = min(RALIGN_CUTOFF,suff_ralign)
      }
    }
  }
  
  ## ----- Primary Criterion 3) Transcriptome Breadth: Fraction of filtered genes detected. -----
  
  if(!is.null(gene_filter)){  
  
    breadth = colMeans(expr > 0)
  
    BREADTH_CUTOFF = hard_breadth
  
    if (!is.null(zcut)){
    
      BREADTH_CUTOFF = max(mean(breadth) - zcut*sd(breadth), BREADTH_CUTOFF) # SD Sub-Criterion
      BREADTH_CUTOFF = max(median(breadth) - zcut*mad(breadth), BREADTH_CUTOFF) # MAD Sub-Criterion
    
      if(mixture){ # Mixture SD Sub-Criterion
      
        is.multimodal = dip.test(breadth)$p.value < dip_thresh
      
        if(is.multimodal){
          mixmdl = normalmixEM(breadth,k=2)
          component = which(mixmdl$mu %in% max(mixmdl$mu))
          BREADTH_CUTOFF = max(mixmdl$mu[component] - zcut*mixmdl$sigma[component], BREADTH_CUTOFF)
        }
      }
    
      if(!is.null(suff_breadth)){
        BREADTH_CUTOFF = min(BREADTH_CUTOFF,suff_breadth)
      }
    }
  }

  ## ----- Primary Criterion 4) FNR AUC. -----
  
  if(!is.null(pos_controls)){  
    
    # Compute FNR AUC  
    ref.glms = simple_FNR_params(expr = expr, pos_controls = pos_controls)
    AUC = NULL
    for (si in 1:dim(sc.eSet)[2]){
      if(!any(is.na(ref.glms[[si]]))){
        AUC[si] = log(exp(ref.glms[[si]][1] + ref.glms[[si]][2] * AUC_range[2]) + 1)/ref.glms[[si]][2] - log(exp(ref.glms[[si]][1] + ref.glms[[si]][2] * AUC_range[1]) + 1)/ref.glms[[si]][2]
      } else {
        stop("glm fit did not converge")
      }
    }
    
    AUC_CUTOFF = hard_AUC
    
    if (!is.null(zcut)){
      
      AUC_CUTOFF = min(mean(AUC) + zcut*sd(AUC), AUC_CUTOFF) # SD Sub-Criterion
      AUC_CUTOFF = min(median(AUC) + zcut*mad(AUC), AUC_CUTOFF) # MAD Sub-Criterion
      
      if(mixture){ # Mixture SD Sub-Criterion
        
        is.multimodal = dip.test(AUC)$p.value < dip_thresh
        
        if(is.multimodal){
          mixmdl = normalmixEM(AUC,k=2)
          component = which(mixmdl$mu %in% min(mixmdl$mu))
          AUC_CUTOFF = min(mixmdl$mu[component] + zcut*mixmdl$sigma[component], AUC_CUTOFF)
        }
      }
      
      if(!is.null(suff_AUC)){
        AUC_CUTOFF = min(AUC_CUTOFF,suff_AUC)
      }
    }
  }
  
  filtered_nreads = logr < LOGR_CUTOFF
  filtered_ralign = ralign < RALIGN_CUTOFF
  filtered_breadth = breadth < BREADTH_CUTOFF
  filtered_fnr = AUC > AUC_CUTOFF
  
  return(list(filtered_nreads = filtered_nreads,
              filtered_ralign = filtered_ralign,
              filtered_breadth = filtered_breadth,
              filtered_fnr = filtered_fnr))
}

