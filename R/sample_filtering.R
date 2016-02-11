#' Fit Logistic Regression Model of FNR against set of positive control (ubiquitously expressed) genes
#' 
#' @details logit(Probability of False Negative) ~ a + b*(mean log10p1 expression) .
#'  
#' @param expr matrix The data matrix in transcript-proportional units (genes in rows, cells in columns).
#' @param pos_controls A boolean vector indexing positive control genes that will be used to compute false-negative rate characteristics.
#' @param fn_tresh Inclusive threshold for negative detection. Default 0.01.
#' 
#' @return A list of logistic regression coefficients corresponding to glm fits in each sample. If a fit did not converge, the result reported is NA.
#'
simple_FNR_params = function(expr, pos_controls, fn_tresh = 0.01){
  
  # Mean log10p1 expression
  mu_obs = rowMeans(log10(expr[pos_controls,]+1))
  
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
#' @param expr matrix The data matrix (genes in rows, cells in columns).
#' @param nreads A numeric vector representing number of reads in each library.
#' If NULL, filtered_nreads will be returned NA.
#' @param ralign A numeric vector representing the proportion of reads aligned to the reference genome in each library.
#' If NULL, filtered_ralign will be returned NA.
#' @param gene_filter A boolean vector indexing genes that will be used to compute library transcriptome breadth.
#' If NULL, filtered_breadth will be returned NA.
#' @param pos_controls A boolean vector indexing positive control genes that will be used to compute false-negative rate characteristics.
#' If NULL, filtered_fnr will be returned NA.
#' @param scale. logical. Will expression be scaled by total expression for FNR computation? Default = FALSE 
#' @param glen Gene lengths for gene-length normalization (normalized data used in FNR computation).
#' @param AUC_range An array of two values, representing range over which FNR AUC will be computed (log10(expr_units + 1)). Default c(0,6)
#' @param zcut A numeric value determining threshold Z-score for sd, mad, and mixture sub-criteria. Default 1.
#' If NULL, only hard threshold sub-criteria will be applied.
#' @param mixture A boolean value determining whether mixture modeling sub-criterion will be applied per primary criterion (metric).
#' If true, a dip test will be applied to each metric. If a metric is multimodal, it is fit to a two-component nomal mixture model. 
#' Samples deviating zcut sd's from optimal mean (in the inferior direction), have failed this sub-criterion.
#' @param dip_thresh A numeric value determining dip test p-value threshold. Default 0.05.
#' @param hard_nreads numeric. Hard (lower bound on) nreads threshold. Default 25000.
#' @param hard_ralign numeric. Hard (lower bound on) ralign threshold. Default 15.
#' @param hard_breadth numeric. Hard (lower bound on) breadth threshold. Default 0.2.
#' @param hard_fnr numeric. Hard (upper bound on) fnr threshold. Default 3.
#' @param suff_nreads numeric. If not null, serves as an upper bound on nreads threshold.
#' @param suff_ralign numeric. If not null, serves as an upper bound on ralign threshold.  Default 65.
#' @param suff_breadth numeric. If not null, serves as an upper bound on breadth threshold.  Default 0.8.
#' @param suff_fnr numeric. If not null, serves as an lower bound on fnr threshold.
#' @param plot_dir If not null, specifies path to plot output
#' @param hist_breaks hist() breaks argument
#' 
#' @return A list with the following elements:
#' \itemize{
#' \item{filtered_nreads}{Boolean. Sample has too few reads.}
#' \item{filtered_ralign}{Boolean. Sample has too few reads aligned.}
#' \item{filtered_breadth}{Boolean. Samples has too few genes detected (low breadth).}
#' \item{filtered_fnr}{Boolean. Sample has a high FNR AUC.}
#' }
#'
#'@importFrom mixtools normalmixEM
#'@importFrom diptest dip.test
#'@export
#'
#'
metric_sample_filter = function(expr, nreads = NULL, ralign = NULL,
                                gene_filter = NULL, pos_controls = NULL,scale. = FALSE,glen = NULL,
                                AUC_range = c(0,6), zcut = 1,
                                mixture = TRUE, dip_thresh = 0.05, 
                                hard_nreads = 25000, hard_ralign = 15, hard_breadth = 0.2, hard_fnr = 3,
                                suff_nreads = NULL, suff_ralign = 65, suff_breadth = 0.8, suff_fnr = NULL,
                                plot_dir = NULL, hist_breaks = 10){
    
  # Create plot directory, if necessary
  if (!is.null(plot_dir) && !file.exists(plot_dir)){
    dir.create(plot_dir)
  }
  
  criterion_count = 0
  
  ### ----- Primary Criterion 1) Number of Reads. -----
  
  if(!is.null(nreads)){ 
    criterion_count = 1
        
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
    filtered_nreads = logr < LOGR_CUTOFF
  }else{
    filtered_nreads = NA
  }
  
  ## ----- Primary Criterion 2) Ratio of reads aligned. -----
  
  if(!is.null(ralign)){
    criterion_count = criterion_count + 1
    
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
    filtered_ralign = ralign < RALIGN_CUTOFF
  }else{
    filtered_ralign = NA
  }
  
  ## ----- Primary Criterion 3) Transcriptome Breadth: Fraction of filtered genes detected. -----
  
  if(!is.null(gene_filter)){  
    criterion_count = criterion_count + 1
  
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
    filtered_breadth = breadth < BREADTH_CUTOFF
  }else{
    filtered_breadth = NA
  }

  ## ----- Primary Criterion 4) FNR AUC. -----
  
  if(!is.null(pos_controls)){  
    criterion_count = criterion_count + 1
    
    # Normalize matrix for FNR estimation
    if(scale.){
      nexpr = mean(colSums(expr))*t(t(expr)/colSums(expr))
    }else{
      nexpr = expr
    }
    if(!is.null(glen)){
      nexpr = mean(glen)*nexpr/glen
    }
        
    # Compute FNR AUC  
    ref.glms = simple_FNR_params(expr = nexpr, pos_controls = pos_controls)
    AUC = NULL
    for (si in 1:dim(expr)[2]){
      if(!any(is.na(ref.glms[[si]]))){
        AUC[si] = log(exp(ref.glms[[si]][1] + ref.glms[[si]][2] * AUC_range[2]) + 1)/ref.glms[[si]][2] - log(exp(ref.glms[[si]][1] + ref.glms[[si]][2] * AUC_range[1]) + 1)/ref.glms[[si]][2]
      } else {
        stop("glm fit did not converge")
      }
    }
    
    AUC_CUTOFF = hard_fnr
    
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
      
      if(!is.null(suff_fnr)){
        AUC_CUTOFF = max(AUC_CUTOFF,suff_fnr)
      }
    }
    filtered_fnr = AUC > AUC_CUTOFF
  }else{
    filtered_fnr = NA
  }
    
  ## ----- Plotting -----
  
  if (!is.null(plot_dir)){
    
      
      pdf(paste0(plot_dir,"/filtering_per_criterion.pdf"))
      
      is_bad = rep(FALSE,dim(expr)[2])
      
      par(mfcol = c(criterion_count,2))
    
      if(!is.null(nreads)){
        is_bad = !filtered_nreads 
        hist(logr, main = paste0("nreads: Thresh = ",signif(LOGR_CUTOFF,3)," , Rm = ",sum(!filtered_nreads)), xlab = "log10(NREADS+1)", breaks = hist_breaks)
        abline(v = LOGR_CUTOFF, col = "red", lty = 2)
      }
    
      if(!is.null(ralign)){
        is_bad = is_bad | !filtered_ralign
        hist(ralign, main = paste0("ralign: Thresh = ",signif(RALIGN_CUTOFF,3)," , Rm = ",sum(!filtered_ralign)), xlab = "RALIGN", breaks = hist_breaks)
        abline(v = RALIGN_CUTOFF, col = "red", lty = 2)
      }
    
      if(!is.null(gene_filter)){
        is_bad = is_bad | !filtered_breadth
        hist(breadth, main = paste0("breadth: Thresh = ",signif(BREADTH_CUTOFF,3)," , Rm = ",sum(!filtered_breadth)), xlab = "BREADTH", breaks = hist_breaks)
        abline(v = BREADTH_CUTOFF, col = "red", lty = 2)
      }
      
      if(!is.null(pos_controls)){
        is_bad = is_bad | !filtered_fnr
        hist(AUC, main = paste0("auc: Thresh = ",signif(AUC_CUTOFF,3)," , Rm = ",sum(!filtered_fnr)), xlab = "FNR AUC", breaks = hist_breaks)
        abline(v = AUC_CUTOFF, col = "red", lty = 2)
      }
      
      if(!is.null(nreads)){
        hist(logr[!is_bad], main = paste0("nreads: Thresh = ",signif(LOGR_CUTOFF,3)," , Rm = ",sum(!filtered_nreads)), xlab = "log10(NREADS+1)", breaks = hist_breaks)
      }
      if(!is.null(ralign)){
        hist(ralign[!is_bad], main = paste0("ralign: Thresh = ",signif(RALIGN_CUTOFF,3)," , Rm, = ",sum(!filtered_ralign)), xlab = "RALIGN", breaks = hist_breaks)
      }
      if(!is.null(gene_filter)){
        hist(breadth[!is_bad],  main = paste0("breadth: Thresh = ",signif(BREADTH_CUTOFF,3)," , Rm = ",sum(!filtered_breadth)), xlab = "BREADTH", breaks = hist_breaks)
      }
      if(!is.null(pos_controls)){
        hist(AUC[!is_bad],  main = paste0("auc: Thresh = ",signif(AUC_CUTOFF,3)," , Rm = ",sum(!filtered_fnr)," , Tot_Rm = ",sum(is_bad)), xlab = "FNR AUC", breaks = hist_breaks)
      }
      dev.off()
    
      pdf(paste0(plot_dir,"/overlap_of_criteria.pdf"))
      v = rbind(!filtered_nreads,!filtered_ralign,!filtered_breadth,!filtered_fnr)
      rownames(v) = c("nreads","ralign","breadth","fnr")
      v = na.omit(v)
      m = v %*% t(v)
      m = t(t(m)/diag(m))
      barplot(m, beside = T,legend.text = T, ylim = c(0,1.5))
      dev.off()
  }
  
  return(list(filtered_nreads = filtered_nreads,
              filtered_ralign = filtered_ralign,
              filtered_breadth = filtered_breadth,
              filtered_fnr = filtered_fnr))
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
#' @param expr matrix The data matrix (genes in rows, cells in columns).
#' @param qual matrix Quality metric data matrix (cells in rows, metrics in columns).
#' @param gene_filter A boolean vector indexing genes that will be used for PCA.
#' If NULL, all genes are used.
#' @param max_exp_pcs numeric number of expression PCs used in quality metric selection. Default 5.
#' @param qual_select_q_thresh numeric. q-value threshold for quality/expression correlation significance tests. Default 0.01
#' @param force_metrics boolean. If not NULL, indexes quality metric to be forcefully included in quality PCA.
#' @param good_metrics boolean. If not NULL, indexes quality metric that indicate better quality when of higher value.
#' @param min_qual_variance numeric. Minimum proportion of selected quality variance addressed in filtering. Default 0.70 
#' @param zcut A numeric value determining threshold Z-score for sd, mad, and mixture sub-criteria. Default 1.
#' @param mixture A boolean value determining whether mixture modeling sub-criterion will be applied per primary criterion (quality score).
#' If true, a dip test will be applied to each quality score. If a metric is multimodal, it is fit to a two-component nomal mixture model. 
#' Samples deviating zcut sd's from optimal mean (in the inferior direction), have failed this sub-criterion.
#' @param dip_thresh A numeric value determining dip test p-value threshold. Default 0.05.
#' @param plot_dir If not null, specifies path to plot output
#' @param hist_breaks hist() breaks argument
#' 
#' @return A boolean, representing samples passing factor-based filter.
#' 
#' @importFrom mixtools normalmixEM
#' @importFrom diptest dip.test
#' @import gplots
#' @export
factor_sample_filter = function(expr, qual, gene_filter = NULL, max_exp_pcs = 5,
                                qual_select_q_thresh = 0.01, force_metrics = NULL, good_metrics = NULL,
                                min_qual_variance = 0.7, zcut = 1, 
                                mixture = TRUE, dip_thresh = .01,
                                plot_dir = NULL, hist_breaks = 20){
    
  # Create plot directory, if necessary
  if (!is.null(plot_dir) && !file.exists(plot_dir)){
    dir.create( plot_dir)
  }

  # Gene filter vector
  if(is.null(gene_filter)){
    gene_filter = rep(TRUE,dim(expr)[1])
  }
    
  ## ----- Expression PCs (based on high-quality genes) -----
  pc = prcomp(t(log1p(expr[gene_filter,])), center = TRUE, scale = TRUE)
  
  ## ----- Quality PCs -----
  cors = cor(pc$x[,1: max_exp_pcs],qual,method = "spearman")
  z = sqrt((dim(pc$x)[1]-3)/1.06)*(1/2)*log((1+cors)/(1-cors))
  p = 2*pnorm(-abs(z))
  to_keep = p.adjust(p,method = "BH") < qual_select_q_thresh
  to_keep_vec = colSums(matrix(to_keep,nrow =  max_exp_pcs)) > 0
  
  if(is.null(good_metrics)){
    good_metrics = rep(FALSE,dim(qual)[2])
  }
  
  # Introduce forced metrics
  if(!is.null(force_metrics)){
    to_keep_vec =  to_keep_vec | force_metrics
  }
  
  if (sum(to_keep_vec) == 0){
    warning("No quality metrics were selected. All samples pass filter.")
    return(rep(TRUE,dim(expr)[2]))
  }
  
  keep_quals = qual[,to_keep_vec]
  
  qpc = prcomp(keep_quals,center = T,scale. = T)
  
  if(!is.null(plot_dir)){
    pdf(paste0(plot_dir,"/qual_csum_var.pdf"))
    csum = cumsum((qpc$sdev^2)/sum(qpc$sdev^2))
    plot(csum, main = "Cumulative Quality PC Variance", ylab = "Fraction of Total Variance")
    abline(h = min_qual_variance, lty = 2, col = "red")
    dev.off()
  }
  num_qual_pcs = which(csum > min_qual_variance)[1]
  
  if(!is.null(plot_dir)){
    for (i in 1:num_qual_pcs){
      pdf(paste0(plot_dir,paste0("/qc_pc",i,".pdf")))
      par(mfrow = c(2,1))
      hist(qpc$x[,i],breaks = hist_breaks, main = paste0("Distribution of Quality PC ",i), xlab = paste0("Qual PC",i))
      barplot(abs(qpc$rotation[,i]),col = c("red","green")[1 + (qpc$rotation[,i] > 0)], cex.names = .25,horiz = T, las=1, main = "Loadings")
      dev.off()
    }
  }
  
  if(!is.null(plot_dir)){
    pdf(paste0(plot_dir,"/qual_corr_heatmap.pdf"))
    heatmap.2(cor(keep_quals),key.title = "",key.xlab = "Spearman Corr.",density.info = 'none',trace = 'none',margins = c(20,20), cexRow = .7, cexCol = .7)
    dev.off()
  }
  
  # Only perform filtering when zcut is not null
  if (!is.null(zcut)){
    
    # Initializing sample removal vector
    to_remove = rep(F,dim(expr)[2])
    
    # Check if "good" metrics have been selected -> if filtering is signed
    is_signed = sum(to_keep_vec & good_metrics) > 0
    
    if(!is_signed){
      warning("No good metrics were selected. Unsigned filtering applied.")
    }
    
    # Loop over quality PCs until we've covered min_qual_variance of the quality variance
    for ( i in 1:num_qual_pcs){
      
      # Quality Score
      qscore = qpc$x[,i] 
      
      # Signed Quality Score
      if (is_signed){
        qscore = qscore * median(sign(qpc$rotation[,i][good_metrics[to_keep_vec]]))
      }
      
      # Simple thresholds
      CUTOFF = median(qscore) - zcut*mad(qscore)
      CUTOFF = max(mean(qscore) - zcut*sd(qscore), CUTOFF)
      
      # Reverse thresholds
      if(!is_signed){
        RCUTOFF = median(qscore) + zcut*mad(qscore)
        RCUTOFF = min(mean(qscore) + zcut*sd(qscore), RCUTOFF)
      }
      
      # Mixture model threshold (Signed Only!)
      if(mixture && is_signed){    
        
        is.multimodal = dip.test(qscore)$p.value < dip_thresh   
        
        if(is.multimodal){
          mixmdl = normalmixEM(qscore,k=2)
          component = which(mixmdl$mu %in% max(mixmdl$mu))
          CUTOFF = max(mixmdl$mu[component] - zcut*mixmdl$sigma[component], CUTOFF)
          
        }
      }
      if(!is.null(plot_dir)){
        
        pdf(paste0(plot_dir,paste0("/qc_filter_pc",i,".pdf")),width = 20)
        par(mfrow = c(1,2))
        hist(qscore,breaks = hist_breaks,main = paste0("Distribution of Quality Score ",i),  xlab = paste0("Qual Score ",i))
        abline(v = CUTOFF, lty = 2, col = "red")
        if(!is_signed){
          abline(v = RCUTOFF, lty = 2, col = "red")
        }
      
        if(i == num_qual_pcs){ # Exceeds minimum quality of variance - last filter
          
          if(is_signed){
            to_remove = to_remove | (qscore < CUTOFF)
            hist(qpc$x[,i][!(qscore < CUTOFF)],breaks = hist_breaks, main = paste0("Thresh = ",signif(CUTOFF,3),", Rm = ",sum(qscore < CUTOFF),", Tot Rm = ",sum(to_remove) ),xlab = paste0("Qual Score ",i))
          }else{
            to_remove = to_remove | ((qscore < CUTOFF) | (qscore > RCUTOFF))
            hist(qpc$x[,i][!((qscore < CUTOFF) | (qscore > RCUTOFF))],breaks = hist_breaks, main = paste0("Threshs = ",signif(CUTOFF,3),"/",signif(RCUTOFF,3),", Rm = ",sum((qscore < CUTOFF) | (qscore > RCUTOFF)),", Tot Rm = ",sum(to_remove) ),xlab = paste0("Qual Score ",i))
          }
          
          dev.off()
          break()
          
        }else{
          if(is_signed){
            to_remove = to_remove | (qscore < CUTOFF)
            hist(qpc$x[,i][!(qscore < CUTOFF)],breaks = hist_breaks, main = paste0("Thresh = ",signif(CUTOFF,3),", Rm = ",sum(qscore < CUTOFF) ),xlab = paste0("Qual Score ",i))
          }else{
            to_remove = to_remove | ((qscore < CUTOFF) | (qscore > RCUTOFF))
            hist(qpc$x[,i][!((qscore < CUTOFF) | (qscore > RCUTOFF))],breaks = hist_breaks, main = paste0("Threshs = ",signif(CUTOFF,3),"/",signif(RCUTOFF,3),", Rm = ",sum((qscore < CUTOFF) | (qscore > RCUTOFF)) ),xlab = paste0("Qual Score ",i))
          }
          dev.off()
        }
      }
      
    }
    
    return(!to_remove)
  }else{
    warning("No z-cutoff specified, thus no filtering results returned.")
    return(NA)
  }
}

