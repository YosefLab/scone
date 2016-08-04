#' Wrapper for Running Essential SCONE Modules
#'
#' @param expr matrix. The expression data matrix (genes in rows, cells in columns).
#' @param qc data frame. The quality control (QC) matrix (cells in rows, metrics in columns)
#'   to be used for filtering, normalization, and evaluation.
#' @param bio factor. The biological condition to be modeled in the Adjustment Step
#'   as variation to be preserved. If adjust_bio="no", it will not be used
#'   for normalization, but only for evaluation.
#' @param batch factor. The known batch variable to be included in the
#'   adjustment model as variation to be removed. If adjust_batch="no", it will
#'   not be used for normalization, but only for evaluation.
#' @param negcon character. The genes to be used as negative controls for filtering, 
#'   normalization, and evaluation. These genes should be expressed uniformily across
#'   the biological phenomenon of interest.
#'   Default NULL.
#' @param verbose character. Verbosity level: higher level is more verbose.
#'   Default "0".
#' @param out_dir character. Output directory. 
#'   Default getwd().
#' @param seed numeric. Random seed. 
#'   Default 112233.
#'
#' @param filt_cells logical. Should cells be filtered? Set to FALSE if low quality cells have already been excluded. If cells are not filtered, then initial gene filtering (the one that is done prior to cell filtering) is disabled as it becomes redundant with the gene filtering that is done after cell filtering.
#'   Default TRUE.
#'      
#' @param filt_genes logical. Should genes be filtered post-sample filtering?
#'   Default TRUE.
#'   
#' @param fnr_maxiter numeric. Maximum number of iterations in EM estimation of expression posteriors.
#'    Dafault 1000.
#'  
#' @param norm_impute character. Should imputation be included in the comparison? 
#'   If 'force', only imputed normalizations will be run.
#'   Default "yes."
#' @param norm_scaling character. Scaling options to be included in the Scaling Step. 
#'   Default c("none", "sum", "deseq", "tmm", "uq", "fq", "detect"). See details.
#' @param norm_rezero logical. Restore imputed zeroes to zero following the Scaling Step?
#'   Default TRUE.
#' @param norm_k_max numeric. Maximum number (norm_k_max) of factors of unwanted variation modeled
#'   in the Adjustment Step.
#'   Default NULL.
#' @param norm_qc_expl numeric. In automatic selection of norm_k_max, what
#'   fraction of variation must be explained by the first norm_k_max PCs of qc? 
#'   Default 0.5. Ignored if norm_k_max is not NULL.
#' @param norm_adjust_bio character. If 'no' it will not be included in the model; if
#'   'yes', both models with and without 'bio' will be run; if 'force', only
#'   models with 'bio' will be run.
#'   Default "yes."
#' @param norm_adjust_batch character. If 'no' it will not be modeled in the Adjustment Step;
#'   if 'yes', both models with and without 'batch' will be run; if 'force',
#'   only models with 'batch' will be run.
#'   Default "yes."
#'   
#' @param eval_dim numeric. The number of principal components to use for evaluation.
#'   Default NULL.
#' @param eval_expr_expl numeric. In automatic selection of eval_dim, what
#'   fraction of variation must be explained by the first eval_dim PCs of expr? 
#'   Default 0.1. Ignored if eval_dim is not NULL.
#' @param eval_poscon character. The genes to be used as positive controls for
#'   evaluation. These genes should be expected to change according to the
#'   biological phenomenon of interest.
#' @param eval_max_kclust numeric. The max number of clusters (> 1) to be used for pam
#'   tightness evaluation. If NULL, tightness will be returned NA.
#' @param eval_stratified_pam logical. If TRUE then maximum ASW for PAM_SIL is
#'   separately computed for each biological-cross-batch condition (accepting
#'   NAs), and a weighted average is returned as PAM_SIL.
#'   Default TRUE.
#' @param report_num numeric. Number of top methods to report.
#'   Default 13.
#'
#' @param out_rda logical. If TRUE, a sconeResults.Rda file with the object that the scone function returns is saved in the out_dir (may be very large for large datasets, but useful for post-processing)
#'   Default FALSE.
#'   
#' @export
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics legend lines points
#' @importFrom utils head sessionInfo write.table
#' @importFrom RColorBrewer brewer.pal
#' @importFrom boot inv.logit
#' @importFrom hexbin plot hexbin
#'
#' @details "ADD DESCRIPTION"
#'
#' @return Directory structure "ADD DESCRIPTION"
#'

scone_easybake <- function(expr, qc,
                  bio=NULL, batch=NULL, negcon=NULL,
                  verbose=c("0","1","2"), out_dir=getwd(), seed = 112233,
                  filt_cells=TRUE,
                  filt_genes=TRUE,
                  
                  fnr_maxiter=1000,
                  
                  norm_impute=c("yes", "no", "force"),
                  norm_scaling=c("none", "sum", "deseq", "tmm", "uq", "fq", "detect"),
                  norm_rezero=TRUE, norm_k_max=NULL, norm_qc_expl=0.5,
                  norm_adjust_bio=c("yes", "no", "force"),
                  norm_adjust_batch=c("yes", "no", "force"), 
                  
                  eval_dim=NULL, eval_expr_expl=0.1,
                  eval_poscon=NULL, eval_max_kclust = 10, eval_stratified_pam = TRUE,
                  
                  report_num = 13, out_rda = FALSE) {

  # Require qc or default (colSums() and colSums(>0))
  verbose = match.arg(verbose)
  verbose = as.numeric(verbose)
  
  printf <- function(...) cat(sprintf(...))
  set.seed(seed)
  
  if(is.null(expr)){stop("expr must be specified")}
  if(is.null(qc)){stop("qc must be specified")}
  if(is.null(negcon)){stop("negcon must be specified")}
  
  if(!is.matrix(expr)){stop("expr must be a matrix")}
  if(!is.data.frame(qc)){stop("qc must be a data frame")}
  
  # Primary Directory
  if (!file.exists(out_dir)) {
    dir.create(out_dir, recursive=TRUE)
  }
  
  # Misc Directory
  misc_dir = paste0(out_dir,"/misc")
  if (!file.exists(misc_dir)) {
    dir.create(misc_dir)
  }
  
  ## ------ Data Filtering Module ------
  if(filt_cells) {
    #Note that there's no need to do the initial gene filtering if
    #there's no cell filtering - it becomes redundant with the final gene filtering
    #(and even makes the final gene filtering stricter than intended because we then do gene filter twice on the same set of cells - increasing the quantile between the first and second iteration)
    
    # Initial Gene Filtering: Select "common" transcripts.
    thresh_fail = quantile(expr[expr > 0], 0.2) ###Make argument###
    num_fail = 10 ###Make argument###
    init.gf.vec = rowSums(expr > thresh_fail) > num_fail
    if(verbose > 0){printf("Data Filtering Module: Initial filter - Kept only %d genes expressed in more than %.2f units in more than %d cells, excluded %d genes\n", sum(init.gf.vec), thresh_fail, num_fail, sum(!init.gf.vec))}
    
    # Initial-Filtering Negative Controls
    small_negcon = intersect(negcon,rownames(expr)[init.gf.vec])
    if(verbose > 0){printf("Data Filtering Module: %d negative control genes to be used for sample filtering\n", length(small_negcon))}
    
    # Metric-based Filtering and Report
    if(verbose > 0){printf("Data Filtering Module: Filtering samples\n")}
    
    pdf(paste0(misc_dir,"/filter_report.pdf"))
    mfilt = metric_sample_filter(expr, mixture = TRUE, plot = TRUE, hist_breaks = 100,
                                 zcut = 4, ###Make argument###
                                 pos_controls = rownames(expr) %in% small_negcon,
                                 gene_filter = init.gf.vec,
                                 nreads = qc$NREADS,ralign = qc$RALIGN) ###Make argument###
    dev.off()
    mfilt = !apply(simplify2array(mfilt[!is.na(mfilt)]),1,any)
    if(verbose > 0){printf("Data Filtering Module: %d samples to be used in downstream analysis, excluded %d samples\n", sum(mfilt),  sum(!mfilt))}
    
    # Implement sample filter
    expr = expr[,mfilt]
    qc = qc[mfilt,]
    if(!is.null(batch)) {
      batch = droplevels(batch[mfilt])
    }
    if(!is.null(bio)) {
      bio = droplevels(bio[mfilt])
    }
  } else {
    if(verbose > 0){printf("Data Filtering Module: Skipping cell filter...\n")}
  }
  
  if(filt_genes){
    thresh_fail = quantile(expr[expr > 0], 0.2) ###Make argument###
    num_fail = 10 ###Make argument###
    final.gf.vec = rowSums(expr > thresh_fail) > num_fail
    if(verbose > 0){printf("Data Filtering Module: Kept only %d genes expressed in more than %.2f units in more than %d cells, excluded %d genes\n", sum(final.gf.vec), thresh_fail, num_fail, sum(!final.gf.vec))}
  
    # Implement gene filter
    expr = expr[final.gf.vec,]
    negcon = negcon[negcon %in% rownames(expr)]
    if(verbose > 0){printf("Data Filtering Module: Kept only %d negative control genes\n", length(negcon))}
    if(!is.null(eval_poscon)){
      eval_poscon = eval_poscon[eval_poscon %in% rownames(expr)]
      if(verbose > 0){printf("Data Filtering Module: Kept only %d positive control genes\n", length(eval_poscon))}
    }
  } else {
    if(verbose > 0){printf("Data Filtering Module: Skipping gene filter...\n")}
  }

  if(verbose > 0){printf("Data Filtering Module: COMPLETE\n\n")}
  ## ------ False-Negative Rate Inference Module ------
  
  if(verbose > 0){printf("False-Negative Rate Inference Module: Estimating posterior expression probability...\n")}
  fnr_out = estimate_ziber(x = expr, 
                           bulk_model = TRUE, ###Make argument###
                           pos_controls = rownames(expr) %in% negcon,
                           verbose = (verbose > 1), 
                           maxiter = fnr_maxiter)
  if(fnr_out$convergence == 1){printf("False-Negative Rate Inference Module: WARNING - EM did not converge. Try increasing fnr_max_iter\n")}
  if(verbose > 0){printf("False-Negative Rate Inference Module: Producing report...\n")}
  
  # FNR Report
  cc <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"),brewer.pal(9,"Set3"))
  
  logistic_coef = simple_FNR_params(expr = expr, pos_controls = rownames(expr) %in% negcon) # "Non-DE" are expected to be expressed in all cells
  
  pdf(paste0(misc_dir,"/fnr_report.pdf"),width = 12,height = 6)
  par(mfrow=c(1,2))
  
  plot(NULL, main = "False Negative Rate Curves",ylim = c(0,1),xlim = c(0,15), ylab = "Failure Probability", xlab = "Median log-Expression")
  x = (0:150)/10
  AUC = NULL
  for(si in 1:ncol(expr)){
    y = inv.logit(logistic_coef[si,1] + logistic_coef[si,2] * x)
    AUC[si] = sum(y)/10
    if(!is.null(batch)){
      colchoice = cc[batch][si]
    }else{
      colchoice = "black"
    }
    lines(x, y, type = 'l', lwd = 2, col = colchoice)
  }
  
  # Barplot of FNR AUC
  if(!is.null(batch)){
    oAUC = order(AUC)
    sAUC = AUC[oAUC]
    sbatch = batch[oAUC]
    barplot(sAUC[order(sbatch)], col=cc[sbatch][order(sbatch)], border=cc[sbatch][order(sbatch)], main="FNR AUC, colored by batch", ylim = c(0,2*max(AUC)))
    legend("topleft", legend=levels(sbatch), fill=cc, cex=0.4)
  }else{
    barplot(sort(AUC), col="black", border="black", main="FNR AUC",ylim = c(0,2*max(AUC)))
  }
  
  hexbin::plot(hexbin(rowMeans(expr > 0)[!is.na(fnr_out$Beta)],inv.logit(fnr_out$Beta[!is.na(fnr_out$Beta)])), xlab = "Detection Rate", ylab = "Expression Rate")
  hexbin::plot(hexbin(fnr_out$Alpha[1,][!is.na(fnr_out$Beta)],inv.logit(fnr_out$Beta[!is.na(fnr_out$Beta)])), ylab = "Expression Rate", xlab = "Median log-Expression")
  
  if(verbose > 0){printf("False-Negative Rate Inference Module: COMPLETE\n\n")}
  
  dev.off()
  
  ## ------ SCONE Prime Module ------
  
  qcmat = as.matrix(qc)
  
  #Allow non-serial (can register serial from outside for debug)
  #BiocParallel::register(BiocParallel::SerialParam()) ###Make argument###
  
  # Imputation arguments
  norm_impute = match.arg(norm_impute)
  imputation = switch(norm_impute, 
                      yes=list(none=impute_null,expect=impute_expectation),
                      no=list(none=impute_null),
                      force=list(expect=impute_expectation))
  impute_args = list(p_nodrop = fnr_out$p_nodrop, mu = exp(fnr_out$Alpha[1,]))
  
  # Scaling arguments
  norm_scaling = match.arg(norm_scaling,several.ok = TRUE)
  
  SUM_FN = function (ei) 
  {
    sums = colSums(ei)
    eo = t(t(ei)*mean(sums)/sums)
    return(eo)
  }
  
  EFF_FN = function (ei) 
  {
    sums = colSums(ei > 0)
    eo = t(t(ei)*sums/mean(sums))
    return(eo)
  }
  
  scaling=list(none=identity,
               sum = SUM_FN,
               deseq=DESEQ_FN,
               tmm = TMM_FN,
               uq = UQ_FN,
               fq = FQT_FN,
               detect = EFF_FN)[norm_scaling]
  
  # Adjustment arguments
  if(is.null(norm_k_max)){
    qc_sd = prcomp(qcmat,center = TRUE,scale = TRUE)$sd
    norm_k_max = which(cumsum(qc_sd^2)/sum(qc_sd^2) > norm_qc_expl)[1]
  }
  if(verbose > 0){printf("Normalization Module: Adjusting for up to %d factors of unwanted variation\n", norm_k_max)}
  
  # Generate Params (RUN = FALSE)
  if(verbose > 0){printf("Normalization Module: Selecting params:\n", norm_k_max)}
  
  params <- scone(expr, imputation = imputation, impute_args = impute_args,
                  scaling=scaling,
                  k_qc=norm_k_max, k_ruv = norm_k_max,
                  qc=qcmat, ruv_negcon = negcon,
                  adjust_bio = match.arg(norm_adjust_bio), bio = bio,
                  adjust_batch = match.arg(norm_adjust_batch), batch = batch,
                  run=FALSE)
  
  is_screened = ((params$imputation_method == "expect") & (params$scaling_method %in% c("detect"))) |
    ((params$adjust_biology == "bio") & (params$adjust_batch != "batch"))
  
  if(norm_rezero){
    is_screened = is_screened | ((params$imputation_method == "expect") & (params$scaling_method %in% c("none")))
  }
  
  params = params[!is_screened,]
  
  if(verbose > 0){print(params)}
  
  # Evaluation arguments
  if(is.null(eval_dim)){
    expr_sd = prcomp(t(expr),center = TRUE,scale = TRUE)$sd
    eval_dim = which(cumsum(expr_sd^2)/sum(expr_sd^2) > eval_expr_expl)[1]
  }
  if(verbose > 0){printf("Normalization Module: Evaluation based on first %d PCs of expression\n", eval_dim)}

  if(is.null(eval_max_kclust)){
    eval_kclust = NULL
    if(verbose > 0){printf("Normalization Module: No PAM-based evaluation\n")}
  }else{
    if(eval_stratified_pam){
      if(!is.null(bio) && !is.null(batch)) {
        temp_lim = min(table(bio,batch))
      } else if(is.null(batch) && is.null(bio)) {
        temp_lim = ncol(expr)
      } else if (!is.null(bio)) { #one of them is not null and the other is
        temp_lim = min(table(bio))
      } else {
        temp_lim = min(table(batch))
      }
      uplim = min(temp_lim-1,eval_max_kclust)
    }else{
      uplim = min(ncol(expr)-1,eval_max_kclust)
    }
    if(uplim >=2){
      eval_kclust = 2:uplim
      if(verbose > 0){printf("Normalization Module: Maximum clustering k for PAM set to %d\n", uplim)}
    }else{
      eval_kclust = NULL
      if(verbose > 0){printf("Normalization Module: No PAM-based evaluation\n")}
    }
    
  }

  # Generate Scores and Ranking
  if(verbose > 0){printf("Normalization Module: Scoring Normalizations...\n")}
  tic = proc.time()
  res <- scone(expr, imputation = imputation, impute_args = impute_args,
                                 scaling=scaling,
                                 k_qc=norm_k_max, k_ruv = norm_k_max,
                                 qc=qcmat, ruv_negcon = negcon,
                                 adjust_bio = match.arg(norm_adjust_bio), bio = bio,
                                 adjust_batch = match.arg(norm_adjust_batch), batch = batch,

                                 run=TRUE, params = params, verbose = (verbose > 1),
                                 eval_poscon = eval_poscon, eval_kclust = eval_kclust, 
                                 stratified_pam = eval_stratified_pam ,rezero = norm_rezero)
  toc = proc.time()
  if(verbose > 0) {
    printf("Normalization Module: Scoring Normalizations Done!...\n")
    printf("time elapsed (seconds):\n")
    print(toc-tic)
  }
  
  
  # Generate SCONE Report
  if(verbose > 0){printf("Normalization Module: Producing main report...\n")}
  pdf(paste0(out_dir,"/scone_report.pdf"),width = 6,height = 6)
  
  pc_obj = prcomp(t(na.omit(t(res$scores[,-ncol(res$scores)]))),center = TRUE, scale = FALSE)
  bpob = biplot_colored(pc_obj,-res$scores[,ncol(res$scores)],expand = .4)
  if(any(rownames(bpob) == "none,none,no_uv,no_bio,no_batch")){
    points(bpob["none,none,no_uv,no_bio,no_batch",1],bpob["none,none,no_uv,no_bio,no_batch",2],lwd = 4, col = "red")
  }
  points(head(as.data.frame(bpob),n = report_num),lwd = 2, col = "blue")
  points(head(as.data.frame(bpob),n = 1),lwd = 2, col = "cyan")
  
  dev.off()
    
  # Recomputing top methods
  params_select = head(res$params,n = report_num)
  if(verbose > 0){printf("Normalization Module: Re-computing top normalizations...\n")}
  tic = proc.time()
  res_select <- scone(expr, imputation = imputation, impute_args = impute_args,
                                 scaling=scaling,
                                 k_qc=norm_k_max, k_ruv = norm_k_max,
                                 qc=qcmat, ruv_negcon = negcon,
                                 adjust_bio = match.arg(norm_adjust_bio), bio = bio,
                                 adjust_batch = match.arg(norm_adjust_batch), batch = batch,
                                 
                                 run=TRUE, return_norm = "in_memory", params = params_select, verbose = (verbose > 1),
                                 eval_poscon = eval_poscon, eval_kclust = eval_kclust, 
                                 stratified_pam = eval_stratified_pam ,rezero = norm_rezero)
  toc = proc.time()
  if(verbose > 0) {
    printf("Normalization Module: Recomputing top methods done!...\n")
    printf("time elapsed (seconds):\n")
    print(toc-tic)
  }
  
  if(out_rda) {
    if(verbose > 0){printf("Normalization Module: Writing sconeResults.Rda file of scone's returned object...\n")}
    save(file=file.path(out_dir, "sconeResults.Rda"), res_select)
  }
  
  if(verbose > 0){
    printf("Normalization Module: (Trumpets): the selected method is... %s!\n", rownames(res_select$scores)[1])
    printf("Normalization Module: Writing top normalization to file...\n")
  }
  
  normle = res_select$normalized_data[[1]]
  write.table(x = normle,row.names = TRUE,col.names = TRUE,quote = FALSE,sep = "\t",eol = "\n",file = paste0(out_dir,"/best_scone.txt"))
  
  if(verbose > 0){printf("Normalization Module: Producing normalization-specific reports...\n")}
  
  write.table(x = res_select$scores, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t", eol = "\n", file = file.path(out_dir,"/normalization_scores.txt"))
  
  for(count in seq_along(rownames(res_select$scores))) {
    nom = rownames(res_select$scores)[count]
    
    if(verbose > 0){printf(paste0("Normalization Module:\t(",count,")\t",nom,"\n"))}
    nom_dir = paste0(misc_dir,"/N",count,"_",gsub(",","_",nom))
    if (!file.exists(nom_dir)) {
      dir.create(nom_dir)
    }
    normle = res_select$normalized_data[[nom]]
    write.table(x = normle,row.names = TRUE,col.names = TRUE,quote = FALSE,sep = "\t",eol = "\n",file = paste0(nom_dir,"/norm.txt"))
    if(is.null(bio) & is.null(batch)){
      pdf(paste0(nom_dir,"/norm_report.pdf"),width = 6,height = 6)
      plot(prcomp(t(normle ))$x, col = "black", pch =16)
      dev.off()
    }else if(is.null(bio) & !is.null(batch)){
      pdf(paste0(nom_dir,"/norm_report.pdf"),width = 6,height = 6)
      plot(prcomp(t(normle ))$x, col = cc[batch], pch =16)
      dev.off()
    }else if(!is.null(bio) & is.null(batch)){
      pdf(paste0(nom_dir,"/norm_report.pdf"),width = 6,height = 6)
      plot(prcomp(t(normle ))$x, col = cc[bio], pch =16)
      dev.off()
    }else{
      pdf(paste0(nom_dir,"/norm_report.pdf"),width = 12,height = 6)
      par(mfrow=c(1,2))
      plot(prcomp(t(normle ))$x, col = cc[batch], pch =16)
      plot(prcomp(t(normle ))$x, col = cc[bio], pch =16)
      dev.off()
    }

  }
  
  if(verbose > 0){printf("Normalization Module: COMPLETE\n")}
  
  ## ------ Bye now ------
  
  if(verbose > 0){
    printf("\n=====sessionInfo=====\n\n")
    sessionInfo()
    }
}