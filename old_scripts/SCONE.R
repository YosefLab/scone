## DR: if I understand correctly "nom" is used only for the name of the file to store the data.
## There should be a better way to do this, using the name of the function (the opposite of get()?).

source("/data/yosef/users/daviderisso/YosefCode/packages/RNASeq/OFBIT/SCONE/estimateFNR.R")
source("/data/yosef/users/daviderisso/YosefCode/packages/RNASeq/OFBIT/SCONE/SCONE_DEFAULTS.R")
library(BiocParallel,quietly = TRUE)
library(cluster,quietly = TRUE)
library(fpc,quietly = TRUE)

runAdjustTask = function(task,evaluate_material,design,nested_model, to_store = FALSE){
  
#???  source("~/YosefCode/packages/RNASeq/OFBIT/SCONE/SCONE.R") # library(SCONE)
  
  result = try({
    
    if(file.exists(paste(evaluate_material$out_dir,"/",task$nom,".Rdata",sep = ""))){
      if(to_store){
        load(paste(evaluate_material$out_dir,"/",task$nom,".Rdata",sep = ""))
        score_out = scone_out$evaluation
      
        # Return values
        if(!is.na(score_out)[1]){
          scores = score_out$scores
        }else{
          scores = NA
        }
  
        print(task$nom)
        return(list(scores = scores))
      }else{
        print(task$nom)
        return()
      }
    
    }else{
    
      norm_lout = ADJUST_FN(task$ei,batch = task$batch,
                            bio = task$condition,
                            uv = task$uv,
                            w = task$w,
                            design = design, 
                            nested_model = nested_model)
      norm_out = exp(norm_lout) - 1
      score_out = scoreMethod(e = norm_out ,condition = evaluate_material$condition, batch = evaluate_material$batch,
                            qual_scores = evaluate_material$qual_scores, HK_scores = evaluate_material$HK_scores, 
                            DE_scores = evaluate_material$DE_scores, CC_scores = evaluate_material$CC_scores,
                            dim_eval = evaluate_material$dim_eval, K_clust = evaluate_material$K_clust, K_NN = evaluate_material$K_NN,
                            fnr_pi = evaluate_material$fnr_pi, nom = task$nom)
    
      scone_out = list(exp_mat = norm_out,evaluation = score_out,method = task$nom)
      save(scone_out,file = paste(evaluate_material$out_dir,"/",task$nom,".Rdata",sep = "")) 
      if(to_store){
      
        # Return values
        if(!is.na(score_out)[1]){
          scores = score_out$scores
        }else{
          scores = NA
        }
    
        print(task$nom)
        return(list(scores = scores))
      }else{
        print(task$nom)
        return()
      }
    
    }
  
  })
  
  return(NA)
  
}

runFactorFreeTask = function(task,evaluate_material,emat){
  
  ## ??? source("~/YosefCode/packages/RNASeq/OFBIT/SCONE/SCONE.R") # library(SCONE)
  
  result = try({
    if(file.exists(paste(evaluate_material$out_dir,"/",task$nom,".Rdata",sep = ""))){
      load(paste(evaluate_material$out_dir,"/",task$nom,".Rdata",sep = ""))
      norm_out = scone_out$exp_mat
      score_out = scone_out$evaluation
    
      # Return values
      if(!is.na(score_out)[1]){
        scores = score_out$scores
      }else{
        scores = NA
      }
    
      print(task$nom)
      return(list(eo = norm_out, scores = scores))
    
      }else{
        norm_out = task$FN(emat)

        score_out = scoreMethod(e = norm_out ,condition = evaluate_material$condition, batch = evaluate_material$batch,
                          qual_scores = evaluate_material$qual_scores, HK_scores = evaluate_material$HK_scores, 
                          DE_scores = evaluate_material$DE_scores, CC_scores = evaluate_material$CC_scores,
                          dim_eval = evaluate_material$dim_eval, K_clust = evaluate_material$K_clust, K_NN = evaluate_material$K_NN,
                          fnr_pi = evaluate_material$fnr_pi, nom = task$nom)
  
      scone_out = list(exp_mat = norm_out,evaluation = score_out,method = task$nom)
      save(scone_out,file = paste(evaluate_material$out_dir,"/",task$nom,".Rdata",sep = "")) 
  
      # Return values
      if(!is.na(score_out)[1]){
        scores = score_out$scores
      }else{
        scores = NA
      }
  
      print(task$nom)
      return(list(eo = norm_out, scores = scores))
    }
  
  })
  
  return(NA)
  
}

scoreMethod = function(e,condition = NULL, batch = NULL, 
                       qual_scores = NULL, HK_scores = NULL,
                       DE_scores = NULL, CC_scores = NULL,
                       dim_eval = NULL, K_clust = NULL, K_NN = NULL,
                       fnr_pi = NULL, nom = NULL){
  if(any(is.na(e) | is.infinite(e) | is.nan(e))){
    warning("NA/Inf/NaN Expression Values. Returning NA.")
    return(NA)
  }
  # Project the data: PCA or wPCA
  if(grepl("^IMPUTE",nom)){
    pc_val = prcomp(t(log(e + 1)),center = TRUE,scale. = TRUE)$x[,1:dim_eval]
  }else{
    pc_val = wPCA(log(e + 1),1 - fnr_pi,nu = dim_eval,filt = TRUE)$x
  }
  
  # Max cor with quality scores
  EXP_QPC_COR = max(cor(pc_val,qual_scores,method = "spearman")^2)
  
  # Max cor with housekeeping scores
  if(!is.null(HK_scores)){
    EXP_HK_COR = max(cor(pc_val,HK_scores,method = "spearman")^2)
  }else{
    EXP_HK_COR = NA
  }
  
  # Max cor with differentially expressed scores
  if(!is.null(DE_scores)){
    EXP_DE_COR = max(cor(pc_val,DE_scores,method = "spearman")^2)
  }else{
    EXP_DE_COR = NA
  }
  
  # Max cor with cell cycle scores
  if(!is.null(CC_scores)){
    EXP_CC_COR = max(cor(pc_val,CC_scores,method = "spearman")^2)
  }else{
    EXP_CC_COR = NA
  }
  
  d = as.matrix(dist(pc_val,method = "euclidean"))
  
  # K-NN Condition
  if(!is.null(condition)){
    m <- sapply(seq_len(nrow(d)), function(s) {
      ds <- d[,s]
      ra <- rank(ds)
      return(mean(condition[(ra <= K_NN+1) & (ra != 1)] == condition[s]))
    })
    KNN_BIO = mean(m)
  }else{
    KNN_BIO = NA
  }
  
  # K-NN Batch
  if(!is.null(batch)){
    m <- sapply(seq_len(nrow(d)), function(s) {
      ds <- d[,s]
      ra <- rank(ds)
      return(mean(batch[(ra <= K_NN+1) & (ra != 1)] == batch[s]))
    })
    KNN_BATCH = mean(m)
  }else{
    KNN_BATCH = NA
  }
  
  # PAM Silh
  if(K_clust > 1){
    pam_object = pam(pc_val,k = K_clust)
    PAM_SIL = pam_object$silinfo$avg.width
    clusters = pam_object$clustering
  }else{
    PAM_SIL = NA
    clusters = NA
  }
  scores = c(EXP_QPC_COR,EXP_HK_COR,EXP_DE_COR,EXP_CC_COR,KNN_BIO,KNN_BATCH,PAM_SIL)
  names(scores) = c("EXP_QPC_COR","EXP_HK_COR","EXP_DE_COR","EXP_CC_COR","KNN_BIO","KNN_BATCH","PAM_SIL")
  return(list(scores = scores, pc_val = pc_val, clusters = clusters))
}

## ===== SCONE: Single-Cell Overview of Normalized Expression data =====
# Intermediate wrapper is missing"
# Decide model for the user - describe reason -"
# One batch one condition - user friendly 
# Check rank of design matrix - if nested - we're good - IFF?
# Positive (evaluate) and negative controls (modeled + evaluate)
# show correlation between q and hk - why we don't try all 
# Default number of dim_UV
# Factor free norms and adjust norms as lists = impute or not, positive data or not, only q , only hk
# function 1, function 2, and k
SCONE = function(e,condition = NULL, batch = NULL, 
                 design = c("factorial","nested"), 
                 nested_model = c("fixed","random"),
                 qual = NULL, is_HK = NULL,
                 is_DE = NULL,is_CC = NULL,
                 dim_UV = 5,
                 dim_eval = 3, K_clust = NULL, K_NN = 10, 
                 out_dir = NULL,
                 K_clust_MAX = 10,K_pseudobatch_MAX = 10,
                 factor_free_only = FALSE, to_store = FALSE, n_workers = 10){
  param = MulticoreParam(workers = n_workers)
  ## ----- Set up output directory ----
  if (file.exists(out_dir)){
    #stop("Designated output directory already exists.")
  } else {
    dir.create(out_dir)
  }
  
  ## ----- Imputation Step ----
  print("Imputation Step...")
  HK_default = FALSE
  if(is.null(is_HK)){
    warning("No control genes have been specified, using all genes for FNR estimation.")
    HK_default = TRUE
    is_HK = rep(TRUE ,nrow(e))
  }
  if(file.exists(paste0(out_dir,"/fnr_out.Rdata"))){
    load(paste0(out_dir,"/fnr_out.Rdata"))
  }else{
    fnr_out = estimateFNR(e, bulk_model = TRUE, is.expressed = is_HK)
    save(fnr_out,file = paste0(out_dir,"/fnr_out.Rdata"))
  }
  fnr_mu = exp(fnr_out$Alpha[1,])
  fnr_pi = (e == 0) * fnr_out$Z 
  imputed_e = e + (fnr_mu * fnr_pi)
  print("Complete.")
  
  ## ----- Generate Unwanted Factors ----
  print("Generating Unwanted Factors for Evaluation...")
  if(file.exists(paste0(out_dir,"/evaluate_material.Rdata"))){
    
    load(paste0(out_dir,"/evaluate_material.Rdata"))
    
  }else{
    
  # Factors of alignment quality metrics
  if(is.null(qual)){
    warning("No quality metrics specified. Total counts and efficiency used.")
    numba = colSums(e == 0)
    beta = colMeans(e == 0)
    qual = cbind(numba, log(beta/(1-beta)))
  }
  
  if(is.null(dim_UV)){
    dim_UV = min(ncol(qual), sum(is_HK))
    warning(paste("Number of unwanted factors not specified. Selecting minimum of number of quality features and number of control genes:",dim_UV))
    if(dim_UV >= ncol(e)){
      stop("Number of unwanted factors >= number of samples!")
    }
  }
  
  qual_scores = prcomp(qual, center = TRUE, scale = TRUE)$x[,1:dim_UV]
  
  if(!HK_default){
    HK_scores = wPCA(log(e[is_HK,] + 1), 1 - fnr_pi[is_HK,], nu = dim_UV, filt = TRUE)$x
  }else{
    HK_scores = NULL
  }
  
  if(!is.null(is_DE)){
    DE_scores = wPCA(log(e[is_DE,] + 1), 1 - fnr_pi[is_DE,], nu = dim_eval, filt = TRUE)$x
  }else{
    DE_scores = NULL
  }
  
  if(!is.null(is_CC)){
    CC_scores = wPCA(log(e[is_CC,] + 1), 1 - fnr_pi[is_CC,], nu = dim_eval, filt = TRUE)$x
  }else{
    CC_scores = NULL
  }
  
  ## ----- Evaluation Criteria ----
  if(is.null(K_clust)){
    if(is.null(DE_scores)){
      warning("No clustering K specified. Selecting K using PAM average silhouette width on raw wPCA.")
      pc_val = wPCA(log(e + 1), 1 - fnr_pi, nu = dim_eval, filt = TRUE)$x
    }else{
      warning("No clustering K specified. Selecting K using PAM average silhouette width on DE scores.")
      pc_val = DE_scores
    }
    
    pamk.best <- pamk(data = pc_val, krange = 1:K_clust_MAX)
    K_clust = pamk.best$nc
  }
  
  data_dir = paste0(out_dir,"/data/")
  if (file.exists(data_dir)){
    #stop("Designated data directory already exists.")
  } else {
    dir.create(data_dir)
  }
  
  evaluate_material = list(condition = condition, batch = batch, 
                           qual_scores = qual_scores, HK_scores = HK_scores,
                           DE_scores = DE_scores, CC_scores = CC_scores,
                           dim_eval = dim_eval, K_clust = K_clust, K_NN = K_NN,
                           fnr_pi = fnr_pi, out_dir = data_dir)
  
  save(evaluate_material,file = paste0(out_dir,"/evaluate_material.Rdata"))
  
  }
  
  print("Complete.")
  
  
  ## ----- Factor-Free Normalization
  ## DR: major change: the current code will copy over and over the same matrices (e and imputed_e).
  ## There's no need to do so as we can infer what matrix to use from the name.
  
  task_list = list()
  
  # Generate Task List for Non-imputed
  task_names = paste("NOIMPUTE",c("NONE","UQ","UQP","DESEQ","DESEQP","FQ","FQP","TMM"),sep = "_")
  
  ##DR: why not TMM__FN_POS?
  task_FN = c(identity,UQ_FN,UQ_FN_POS,DESEQ_FN,DESEQ_FN_POS,FQ_FN,FQ_FN_POS,TMM_FN)
  names(task_FN) = task_names
  for (task in task_names){
    task_list[[task]] = list(FN = task_FN[[task]], nom = task)
  }
  # Generate Task List for Imputed
  task_names = paste("IMPUTE",c("NONE","UQ","DESEQ","FQ","TMM"),sep = "_")
  task_FN = c(identity,UQ_FN,DESEQ_FN,FQ_FN,TMM_FN)
  names(task_FN) = task_names
  for (task in task_names){
    task_list[[task]] = list(FN = task_FN[[task]], nom = task)
  }
  
  # Parallel Normalization + Evaluation
  print("Factor-Free Normalization and Evaluation...")
  factor_free_out_noimpute = bplapply(task_list[grep("^NOIMPUTE", names(task_list))], FUN = runFactorFreeTask, evaluate_material = evaluate_material, emat=e, BPPARAM = param)
  factor_free_out_impute = bplapply(task_list[grep("^IMPUTE", names(task_list))], FUN = runFactorFreeTask, evaluate_material = evaluate_material, emat=imputed_e, BPPARAM = param)
  factor_free_out <- c(factor_free_out_noimpute, factor_free_out_impute)
  print("Complete.")
  
  ## ----- Selection Step (Optional)
  if(factor_free_only){
    return(list(factor_free_out = factor_free_out, factor_based_out = NA))
  }
  ## ----- Factor-Based Normalization
  
  print("Selecting Base Tasks...")
  base_task_list <- lapply(1:length(factor_free_out), function(i) {
    if(!any(is.na(factor_free_out[[i]]$eo) | is.nan(factor_free_out[[i]]$eo) | is.infinite(factor_free_out[[i]]$eo))) {
      return(list(ei = factor_free_out[[i]]$eo, nom = names(factor_free_out)[i]))
    }
  })
  names(base_task_list) <- names(factor_free_out)
  base_task_list <- base_task_list[!sapply(base_task_list, is.null)]
  print("Complete.")
  
  
  # WEIGHT
  print("Extending Tasks by Weighting Scheme...")
  ext_task_list = list()
  for( task in names(base_task_list)){
    base_nom = task
    if(!grepl("^IMPUTE",task)){
      ext_nom = paste(base_nom,"ZEROWEIGHT",sep = "_")
      ext_task_list[[ext_nom]] = base_task_list[[base_nom]]
      ext_task_list[[ext_nom]]$nom = ext_nom 
      ext_task_list[[ext_nom]]$w = 0 + (ext_task_list[[ext_nom]]$ei > 0)
    }
    ext_nom = paste(base_nom,"NOWEIGHT",sep = "_")
    ext_task_list[[ext_nom]] = base_task_list[[base_nom]]
    ext_task_list[[ext_nom]]$nom = ext_nom 
    ext_task_list[[ext_nom]]$w = NULL
  }
  print("Complete.")
  
  base_task_list = ext_task_list
  
  # Generate Method-Specific Factors (Should be parallel)
  print("Computing Method-Specific Factors...")
  
  task_factor_list = list()  
  for( task in names(base_task_list)){
    
    # HK
    if(!HK_default){
      if(!is.null(base_task_list[[task]]$w)){
        task_factor_list[[task]]$HK_scores = wPCA(x = log(base_task_list[[task]]$ei[is_HK,]+1),w = base_task_list[[task]]$w[is_HK,],nu = dim_UV,filt = TRUE)$x
      }else{
        task_factor_list[[task]]$HK_scores = prcomp(t(log(base_task_list[[task]]$ei[is_HK,]+1)))$x[,1:dim_UV]
      }
      
      if(is.null(K_pseudobatch)){
        pamk_object = pamk(task_factor_list[[task]]$HK_scores,krange = 1:K_pseudobatch_MAX)
        task_factor_list[[task]]$hk_nc = pamk_object$nc
        task_factor_list[[task]]$pam_clustering =  pamk_object$pamobject$clustering
      }else{
        task_factor_list[[task]]$hk_nc = K_pseudobatch
        task_factor_list[[task]]$pam_clustering = pam(task_factor_list[[task]]$HK_scores,k = K_pseudobatch)$clustering
      }
    }else{
      task_factor_list[[task]]$HK_scores = NULL
      task_factor_list[[task]]$hk_nc = NULL
      task_factor_list[[task]]$pam_clustering = NULL
    }
  }
  print("Complete.")
  
  # Generate General Factor Clusterings
  print("Generating General Factor-Based Clusters...")
  if(is.null(K_pseudobatch)){
    pamk_object = pamk(evaluate_material$qual_scores,krange = 1:K_pseudobatch_MAX)
    qual_nc = pamk_object$nc
    qual_clustering =  as.factor(pamk_object$pamobject$clustering)
  }else{
    qual_nc = K_pseudobatch
    qual_clustering = as.factor(pam(evaluate_material$qual_scores,k = K_pseudobatch)$clustering)
  }
  print("Complete.")
  
  # BIO
  print("Extending Tasks by Biological Covariate...")
  ext_task_list = list()
  for( task in names(base_task_list)){
    base_nom = task
    
    if(!is.null(condition)){
      ext_nom = paste(base_nom,"BIO",sep = "_")
      ext_task_list[[ext_nom]] = base_task_list[[base_nom]]
      ext_task_list[[ext_nom]]$nom = ext_nom 
      ext_task_list[[ext_nom]]$condition = condition
    }
    ext_nom = paste(base_nom,"NOBIO",sep = "_")
    ext_task_list[[ext_nom]] = base_task_list[[base_nom]]
    ext_task_list[[ext_nom]]$nom = ext_nom 
    ext_task_list[[ext_nom]]$condition = NULL
  }
  print("Complete.")
  
  print("Extending Tasks by Batch Covariate...")
  base_task_list = ext_task_list
  # BATCH
  ext_task_list = list()
  for( task in names(base_task_list)){
    base_nom = task
    
    if(!is.null(batch)){
      ext_nom = paste(base_nom,"BATCH",sep = "_")
      ext_task_list[[ext_nom]] = base_task_list[[base_nom]]
      ext_task_list[[ext_nom]]$nom = ext_nom 
      ext_task_list[[ext_nom]]$batch = batch
    }
    ext_nom = paste(base_nom,"NOBATCH",sep = "_")
    ext_task_list[[ext_nom]] = base_task_list[[base_nom]]
    ext_task_list[[ext_nom]]$nom = ext_nom 
    ext_task_list[[ext_nom]]$batch = NULL
  }
  print("Complete.")
  
  print("Extending Tasks by UV Covariates...")
  base_task_list = ext_task_list  
  # UV
  ext_task_list = list()
  for( task in names(base_task_list)){
    base_nom = task
    method_class = gsub("WEIGHT_.*","WEIGHT",task)
    
    ext_nom = paste(base_nom,"NOUV",sep = "_")
    ext_task_list[[ext_nom]] = base_task_list[[base_nom]]
    ext_task_list[[ext_nom]]$nom = ext_nom 
    ext_task_list[[ext_nom]]$uv = NULL
    
    # HK UV
    if(!HK_default){
      for(j in 1:dim_UV){
        ext_nom = paste(base_nom,paste("HK",j,sep = "_"),sep = "_")
        ext_task_list[[ext_nom]] = base_task_list[[base_nom]]
        ext_task_list[[ext_nom]]$nom = ext_nom 
        ext_task_list[[ext_nom]]$uv = as.matrix(task_factor_list[[method_class]]$HK_scores[,1:j])
      }  
      #       if(task_factor_list[[method_class]]$hk_nc > 1){
      #         ext_nom = paste(base_nom,"HKCLUST",sep = "_")
      #         ext_task_list[[ext_nom]] = base_task_list[[base_nom]]
      #         ext_task_list[[ext_nom]]$nom = ext_nom 
      #         ext_task_list[[ext_nom]]$uv = as.matrix(as.factor(task_factor_list[[method_class]]$pam_clustering))
      #       }
    }
    # Qual UV
    for(j in 1:dim_UV){
      ext_nom = paste(base_nom,paste("Q",j,sep = "_"),sep = "_")
      ext_task_list[[ext_nom]] = base_task_list[[base_nom]]
      ext_task_list[[ext_nom]]$nom = ext_nom 
      ext_task_list[[ext_nom]]$uv = as.matrix(evaluate_material$qual_scores[,1:j])
    }  
    #     if(qual_nc > 1){
    #       ext_nom = paste(base_nom,"QCLUST",sep = "_")
    #       ext_task_list[[ext_nom]] = base_task_list[[base_nom]]
    #       ext_task_list[[ext_nom]]$nom = ext_nom 
    #       ext_task_list[[ext_nom]]$uv = as.matrix(as.factor(qual_clustering))
    #     }
  }
  print("Complete.")
  ext_task_list = ext_task_list[!grepl("NOBATCH_NOUV",names(ext_task_list))]
  if(!to_store){
    finished_tasks = gsub(".Rdata","",list.files(evaluate_material$out_dir))
    ext_task_list = ext_task_list[! (names(ext_task_list) %in% finished_tasks)]
  }
  print(length(ext_task_list))
  print("Factor-Based Adjustment Normalization and Evaluation...")
  factor_based_out = bplapply(ext_task_list,FUN = runAdjustTask, evaluate_material = evaluate_material,design = design,nested_model = nested_model,BPPARAM = param)
  print("Complete.")
  
  return(list(factor_free_out = factor_free_out, factor_based_out = factor_based_out))
}
