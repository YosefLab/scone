#' Wrapper for Running Essential SCONE Modules
#' 
#' @param expr matrix. The expression data matrix (genes in rows, cells in
#'   columns).
#' @param qc data frame. The quality control (QC) matrix (cells in rows, 
#'   metrics in columns) to be used for filtering, normalization,
#'   and evaluation.
#' @param bio factor. The biological condition to be modeled in the Adjustment
#'   Step as variation to be preserved. If adjust_bio="no", it will not be used
#'   for normalization, but only for evaluation.
#' @param batch factor. The known batch variable to be included in the 
#'   adjustment model as variation to be removed. If adjust_batch="no", it will
#'   not be used for normalization, but only for evaluation.
#' @param negcon character. The genes to be used as negative controls for
#'   filtering, normalization, and evaluation. These genes should be expressed
#'   uniformily across the biological phenomenon of interest. Default NULL.
#' @param verbose character. Verbosity level: higher level is more verbose. 
#'   Default "0".
#' @param out_dir character. Output directory. Default getwd().
#' @param seed numeric. Random seed. Default 112233.
#'   
#' @param filt_cells logical. Should cells be filtered? Set to FALSE if low
#'   quality cells have already been excluded. If cells are not filtered, then
#'   initial gene filtering (the one that is done prior to cell filtering) is
#'   disabled as it becomes redundant with the gene filtering that is done 
#'   after cell filtering. Default TRUE.
#' @param filt_genes logical. Should genes be filtered post-sample filtering? 
#'   Default TRUE.
#' @param always_keep_genes logical. A character vector of gene names that
#'   should never be excluded (e.g., marker genes). Default NULL.
#' @param fnr_maxiter numeric. Maximum number of iterations in EM estimation of
#'   expression posteriors. If 0, then FNR estimation is skipped entirely, and
#'   as a consequence no imputation will be performed, disregarding the value
#'   of the "norm_impute" argument. Default 1000.
#'   
#' @param norm_impute character. Should imputation be included in the
#'   comparison? If 'force', only imputed normalizations will be run. Default
#'   "yes."
#' @param norm_scaling character. Scaling options to be included in the Scaling
#'   Step. Default c("none", "sum", "deseq", "tmm", "uq", "fq", "detect"). See
#'   details.
#' @param norm_rezero logical. Restore imputed zeroes to zero following the
#'   Scaling Step? Default TRUE.
#' @param norm_k_max numeric. Max number (norm_k_max) of factors of unwanted
#'   variation modeled in the Adjustment Step. Default NULL.
#' @param norm_qc_expl numeric. In automatic selection of norm_k_max, what 
#'   fraction of variation must be explained by the first norm_k_max PCs of qc?
#'   Default 0.5. Ignored if norm_k_max is not NULL.
#' @param norm_adjust_bio character. If 'no' it will not be included in the
#'   model; if 'yes', both models with and without 'bio' will be run; if
#'   'force', only models with 'bio' will be run. Default "yes."
#' @param norm_adjust_batch character. If 'no' it will not be modeled in the
#'   Adjustment Step; if 'yes', both models with and without 'batch' will be
#'   run; if 'force', only models with 'batch' will be run. Default "yes."
#'   
#' @param eval_dim numeric. The number of principal components to use for
#'   evaluation. Default NULL.
#' @param eval_expr_expl numeric. In automatic selection of eval_dim, what 
#'   fraction of variation must be explained by the first eval_dim PCs of expr?
#'   Default 0.1. Ignored if eval_dim is not NULL.
#' @param eval_poscon character. The genes to be used as positive controls for 
#'   evaluation. These genes should be expected to change according to the 
#'   biological phenomenon of interest.
#' @param eval_max_kclust numeric. The max number of clusters (> 1) to be used
#'   for pam tightness evaluation. If NULL, tightness will be returned NA.
#' @param eval_stratified_pam logical. If TRUE then maximum ASW for PAM_SIL is 
#'   separately computed for each biological-cross-batch condition (accepting 
#'   NAs), and a weighted average is returned as PAM_SIL. Default TRUE.
#' @param report_num numeric. Number of top methods to report. Default 13.
#'   
#' @param out_rda logical. If TRUE, sconeResults.Rda file with the object that
#'   the scone function returns is saved in the out_dir (may be very large for
#'   large datasets, but useful for post-processing) Default FALSE.
#' @param eval_negcon character. Alternative negative control gene list for
#'   evaluation only.
#' @param ... extra params passed to the metric_sample_filter and scone when
#'   they're called by easybake
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
#' @examples
#' set.seed(101)
#' mat <- matrix(rpois(1000, lambda = 5), ncol=10)
#' colnames(mat) <- paste("X", 1:ncol(mat), sep="")
#' obj <- SconeExperiment(mat)
#' res <- scone(obj, scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
#'            evaluate=TRUE, k_ruv=0, k_qc=0, eval_kclust=2, 
#'            bpparam = BiocParallel::SerialParam())
#' qc = as.matrix(cbind(colSums(mat),colSums(mat > 0)))
#' rownames(qc) = colnames(mat)
#' colnames(qc) = c("NREADS","RALIGN")
#' \dontrun{
#' scone_easybake(mat, qc = as.data.frame(qc), verbose = "2", 
#'    norm_adjust_bio= "no",
#'    norm_adjust_batch= "no", norm_k_max = 0,
#'    fnr_maxiter = 0, filt_cells=FALSE, filt_genes=FALSE,
#'    eval_stratified_pam = FALSE,
#'    out_dir="~/scone_out")
#' }


scone_easybake <- function(expr, qc,
                           bio=NULL, batch=NULL, negcon=NULL,
                           verbose=c("0","1","2"), 
                           out_dir=getwd(), seed = 112233,
                           filt_cells=TRUE,
                           filt_genes=TRUE, always_keep_genes = NULL,
                           
                           fnr_maxiter=1000,
                           
                           norm_impute=c("yes", "no", "force"),
                           norm_scaling=c("none", "sum", "deseq",
                                          "tmm", "uq", "fq", "detect"),
                           norm_rezero=TRUE, norm_k_max=NULL, norm_qc_expl=0.5,
                           norm_adjust_bio=c("yes", "no", "force"),
                           norm_adjust_batch=c("yes", "no", "force"), 
                           
                           eval_dim = NULL, eval_expr_expl = 0.1,
                           eval_poscon = NULL, eval_negcon = negcon,
                           eval_max_kclust = 10, eval_stratified_pam = TRUE,
                           
                           report_num = 13, out_rda = FALSE, ...) {
  
  # Require qc or default (colSums() and colSums(>0))
  verbose = match.arg(verbose)
  verbose = as.numeric(verbose)
  
  printf <- function(...) cat(sprintf(...))
  set.seed(seed)
  
  if(is.null(expr)){stop("expr must be specified")}
  if(is.null(qc)){stop("qc must be specified")}
  if(is.null(negcon) && fnr_maxiter > 0){
    stop("negcon must be specified if fnr_maxiter > 0")
  }
  
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
  
  #duplicate stdout to log file
  logFile = file(file.path(out_dir, "stdout.txt"), open = "wt")
  sink(logFile, type = "output", split = TRUE)
  
  ## ------ Data Filtering Module ------
  if(filt_cells) {
    
    # Note: There's no need to do the initial gene-filtering if
    # there's no cell filtering - it becomes redundant with the
    # final gene filtering (and even makes the final gene filtering
    # stricter than intended because we then do gene filter twice 
    # on the same set of cells - increasing the quantile between 
    # the first and second iteration)
    
    # Initial Gene Filtering: Select "common" transcripts.
    thresh_fail = quantile(expr[expr > 0], 0.2) ###Make argument###
    num_fail = 10 ###Make argument###
    init.gf.vec = rowSums(expr > thresh_fail) > num_fail
    if(verbose > 0){printf(paste0("Data Filtering Module: Initial filter",
                                  " - Kept only %d genes expressed in more ",
                                  "than %.2f units in more than %d cells, ",
                                  "excluded %d genes\n"), sum(init.gf.vec), 
                           thresh_fail, num_fail, sum(!init.gf.vec))}
    
    # Initial-Filtering Negative Controls
    small_negcon = intersect(negcon,rownames(expr)[init.gf.vec])
    if(verbose > 0){printf(paste0("Data Filtering Module: %d negative ",
                                  "control genes to be used for sample",
                                  " filtering\n"), length(small_negcon))}
    
    # Metric-based Filtering and Report
    if(verbose > 0){printf("Data Filtering Module: Filtering samples\n")}
    
    pdf(paste0(misc_dir,"/filter_report.pdf"))
    mfilt = metric_sample_filter(expr, plot = TRUE, hist_breaks = 100,
                                 zcut = 4, ###Make argument###
                                 pos_controls = rownames(expr) %in%
                                   small_negcon,
                                 gene_filter = init.gf.vec,
                                 nreads = qc$NREADS, 
                                 ralign = qc$RALIGN, ...) ###Make argument###
    dev.off()
    mfilt = !apply(simplify2array(mfilt[!is.na(mfilt)]),1,any)
    if(verbose > 0){printf(paste0("Data Filtering Module: %d samples to be ",
                                  "used in downstream analysis, excluded %d ",
                                  "samples\n"), sum(mfilt),  sum(!mfilt))}
    
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
    if(!is.null(always_keep_genes)) {
      final.gf.vec[always_keep_genes %in% rownames(expr)] = TRUE  
    }
    if(verbose > 0){printf(paste0("Data Filtering Module: Kept only %d genes ",
                                  "expressed in more than %.2f units in more",
                                  " than %d cells (or ones that were forced ",
                                  "to be kept by the user's argument), ",
                                  "excluded %d genes\n"), sum(final.gf.vec), 
                           thresh_fail, num_fail, sum(!final.gf.vec))}
    
    # Implement gene filter
    expr = expr[final.gf.vec,]
    negcon = negcon[negcon %in% rownames(expr)]
    eval_negcon = eval_negcon[eval_negcon %in% rownames(expr)]
    if(verbose > 0){printf(paste0("Data Filtering Module: Kept only %d ",
                                  "negative control genes\n"), length(negcon))}
    if(!is.null(eval_poscon)){
      eval_poscon = eval_poscon[eval_poscon %in% rownames(expr)]
      if(verbose > 0){printf(paste0("Data Filtering Module: ",
                                    "Kept only %d positive control genes\n"),
                             length(eval_poscon))}
    }
  } else {
    if(verbose > 0){printf("Data Filtering Module: Skipping gene filter...\n")}
  }
  
  if(verbose > 0){printf("Data Filtering Module: COMPLETE\n\n")}
  ## ------ False-Negative Rate Inference Module ------
  
  cc <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"),brewer.pal(9,"Set3"))
  
  if(fnr_maxiter > 0) {
    ## fnr_maxiter > 0 -- user wants to do fnr estimation
    
    if(verbose > 0){printf(paste0("False-Negative Rate Inference Module:",
                                  " Estimating posterior ",
                                  "expression probability...\n"))}
    fnr_out = estimate_ziber(x = expr, 
                             bulk_model = TRUE, ###Make argument###
                             pos_controls = rownames(expr) %in%
                               negcon,
                             verbose = (verbose > 1), 
                             maxiter = fnr_maxiter)
    if(fnr_out$convergence == 1){printf(paste0("False-Negative Rate ",
                                               "Inference Module: WARNING - ",
                                               "EM did not converge. Try ",
                                               "increasing fnr_max_iter\n"))}
    if(verbose > 0){printf(paste0("False-Negative Rate ",
                                  "Inference Module: Producing report...\n"))}
    
    # FNR Report
    
    logistic_coef = simple_FNR_params(expr = expr, 
                                      pos_controls = rownames(expr) %in%
                                        negcon)
    # "Non-DE" are expected to be expressed in all cells
    
    pdf(paste0(misc_dir,"/fnr_report.pdf"),width = 12,height = 6)
    par(mfrow=c(1,2))
    
    plot(NULL, main = "False Negative Rate Curves",ylim = c(0,1),
         xlim = c(0,15), ylab = "Failure Probability",
         xlab = "Median log-Expression")
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
      barplot(sAUC[order(sbatch)], col=cc[sbatch][order(sbatch)],
              border=cc[sbatch][order(sbatch)],
              main="FNR AUC, colored by batch",
              ylim = c(0,2*max(AUC)))
      if(fnr_out$convergence != 1) {
        # this line throws an exception 
        # if not converged: "'legend' is of length 0"
        legend("topleft", legend=levels(sbatch), fill=cc, cex=0.4)
      }
    }else{
      barplot(sort(AUC), col="black", border="black",
              main="FNR AUC",ylim = c(0,2*max(AUC)))
    }
    
    hexbin::plot(hexbin(rowMeans(expr > 0)[!is.na(fnr_out$Beta)],
                        inv.logit(fnr_out$Beta[!is.na(fnr_out$Beta)])), 
                 xlab = "Detection Rate", ylab = "Expression Rate")
    hexbin::plot(hexbin(fnr_out$Alpha[1,][!is.na(fnr_out$Beta)],
                        inv.logit(fnr_out$Beta[!is.na(fnr_out$Beta)])),
                 ylab = "Expression Rate", xlab = "Median log-Expression")
    
    if(verbose > 0){printf(paste0("False-Negative Rate ",
                                  "Inference Module: COMPLETE\n\n"))}
    
    dev.off()
    
  } else {
    ## fnr_maxiter <= 0 (should be 0, check for <= 0 for robustness)
    ## -- user doesn't want to do fnr estimation
    if(verbose > 0) {
      printf("Skipping False-Negative Rate Inference Module...\n")
      printf("(and consequently no imputation will be performed too)\n\n")
      norm_impute = "no"
    }
    
  }
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
  if((norm_impute %in% c("yes", "force")) && fnr_maxiter > 0) {
    impute_args = list(p_nodrop = fnr_out$p_nodrop,
                       mu = exp(fnr_out$Alpha[1,]))
  } else {
    impute_args = NULL 
  }
  
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
  if(verbose > 0){printf(paste0("Normalization Module: Adjusting ",
                                "for up to %d factors",
                                " of unwanted variation\n"),
                         norm_k_max)}
  
  # Generate Params (RUN = FALSE)
  if(verbose > 0){printf("Normalization Module: Selecting params:\n")}
  
  args_list = list(object = expr,
                   qc=qcmat)
  
  # Creating a SconeExperiment Object
  if(!is.null(negcon)){ 
    args_list = c( args_list, 
                   list( negcon_ruv = rownames(expr) %in% negcon ))
  }
  if(!is.null( eval_negcon)){ 
    args_list = c( args_list, 
                   list( negcon_eval = rownames(expr) %in% eval_negcon ))
  }
  if(!is.null(eval_poscon)){ 
    args_list = c( args_list, 
                   list( poscon = rownames(expr)%in% eval_poscon ))
  }
  if(!is.null(bio)){
    args_list = c( args_list, list( bio = bio ))
  }
  if(!is.null(batch)){
    args_list = c( args_list, list( batch = batch ))
  }
  
  my_scone <- do.call(SconeExperiment,args_list)
  
  my_scone <- scone(my_scone,
                    imputation = imputation, impute_args = impute_args,
                    scaling=scaling,
                    k_qc=norm_k_max, k_ruv = norm_k_max,
                    
                    adjust_bio = match.arg(norm_adjust_bio),
                    adjust_batch = match.arg(norm_adjust_batch),
                    eval_kclust = NULL,
                    run=FALSE, ...)
  
  is_screened = ((get_params(my_scone)$imputation_method == "expect") &
                   (get_params(my_scone)$scaling_method %in% c("detect"))) |
    ((get_params(my_scone)$adjust_biology == "bio") 
     & (get_params(my_scone)$adjust_batch != "batch"))
  
  if(norm_rezero){
    is_screened = is_screened | 
      ((get_params(my_scone)$imputation_method == "expect") &
         (get_params(my_scone)$scaling_method %in%
            c("none")))
  }
  
  my_scone = select_methods(my_scone, 
                            rownames(get_params(my_scone))[!is_screened ])
  
  if(verbose > 0){print(get_params(my_scone))}
  
  # Evaluation arguments
  if(is.null(eval_dim)){
    expr_sd = prcomp(t(expr),center = TRUE,scale = TRUE)$sd
    eval_dim = which(cumsum(expr_sd^2)/sum(expr_sd^2) > eval_expr_expl)[1]
  }
  if(verbose > 0){printf(paste0("Normalization Module: Evaluation ",
                                "based on first %d PCs of expression\n"),
                         eval_dim)}
  
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
      if(verbose > 0){printf(paste0("Normalization Module: ",
                                    " Maximum clustering k for",
                                    " PAM set to %d\n"),
                             uplim)}
    }else{
      eval_kclust = NULL
      if(verbose > 0){printf(paste0("Normalization Module:",
                                    " No PAM-based evaluation\n"))}
    }
    
  }
  
  
  # Generate Scores and Ranking
  if(verbose > 0){printf("Normalization Module: Scoring Normalizations...\n")}
  tic = proc.time()
  
  my_scone <- scone(my_scone,
                    imputation = imputation, impute_args = impute_args,
                    scaling=scaling,
                    k_qc=norm_k_max, k_ruv = norm_k_max,
                    
                    adjust_bio = match.arg(norm_adjust_bio),
                    adjust_batch = match.arg(norm_adjust_batch),
                    run=TRUE, verbose = (verbose > 1),
                    stratified_pam = eval_stratified_pam, 
                    eval_kclust = eval_kclust,
                    rezero = norm_rezero, ...)
  
  toc = proc.time()
  if(verbose > 0) {
    printf("Normalization Module: Scoring Normalizations Done!...\n")
    printf("time elapsed (seconds):\n")
    print(toc-tic)
  }
  
  
  # Generate SCONE Report
  if(verbose > 0){printf("Normalization Module: Producing main report...\n")}
  pdf(paste0(out_dir,"/scone_report.pdf"),width = 6,height = 6)
  
  pc_obj = prcomp(apply(t(get_scores(my_scone)),1,rank),
                  center = TRUE,scale = FALSE)
  bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)
  
  points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1)
  points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1.5)
  
  points(t(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",]),
         pch = 1, col = "blue", cex = 1)
  points(t(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",]),
         pch = 1, col = "blue", cex = 1.5)
  
  arrows(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",][1],
         bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",][2],
         bp_obj[1,][1],
         bp_obj[1,][2],
         lty = 2, lwd = 2)
  
  dev.off()
  
  # Recomputing top methods
  params_select = head(get_params(my_scone),n = report_num)
  if(verbose > 0){printf(paste0("Normalization Module: ",
                                "Re-computing top normalizations...\n"))}
  tic = proc.time()
  
  scone_res = list()
  
  # List of normalized matrices
  scone_res$normalized_data = lapply(as.list(rownames(params_select)), 
                                     FUN = function(z){
                                       get_normalized(my_scone,
                                                      method = z,log=TRUE)
                                     })
  names(scone_res$normalized_data) = rownames(params_select)
  
  # Parameter matrix
  scone_res$params = get_params(my_scone)[rownames(params_select),]
  
  # Merged score matrix
  scone_res$scores = cbind(get_scores(my_scone),
                           get_score_ranks(my_scone))[rownames(params_select),]
  colnames(scone_res$scores)[ncol(scone_res$scores)] = "mean_score_rank"
  
  
  toc = proc.time()
  if(verbose > 0) {
    printf("Normalization Module: Recomputing top methods done!...\n")
    printf("time elapsed (seconds):\n")
    print(toc-tic)
  }
  
  if(out_rda) {
    if(verbose > 0){printf(paste0("Normalization Module: ",
                                  "Writing sconeResults.Rda file ",
                                  "of scone's returned object...\n"))}
    save(file=file.path(out_dir, "sconeResults.Rda"),  scone_res)
  }
  
  if(verbose > 0){
    printf("Normalization Module: (Trumpets): the selected method is... %s!\n",
           rownames( scone_res$scores)[1])
    printf("Normalization Module: Writing top normalization to file...\n")
  }
  
  normle =  scone_res$normalized_data[[1]]
  write.table(x = normle,row.names = TRUE,col.names = TRUE,quote = FALSE,
              sep = "\t",eol = "\n",file = paste0(out_dir,"/best_scone.txt"))
  
  if(verbose > 0){printf(paste0("Normalization Module: ",
                                "Producing normalization-specific ",
                                "reports...\n"))}
  
  write.table(x =  scone_res$scores, row.names = TRUE, col.names = TRUE,
              quote = FALSE, sep = "\t", eol = "\n",
              file = file.path(out_dir,"/normalization_scores.txt"))
  
  for(count in seq_along(rownames( scone_res$scores))) {
    nom = rownames( scone_res$scores)[count]
    
    if(verbose > 0){
      printf(paste0("Normalization Module:\t(",count,")\t",nom,"\n"))
    }
    nom_dir = paste0(misc_dir,"/N",count,"_",gsub(",","_",nom))
    if (!file.exists(nom_dir)) {
      dir.create(nom_dir)
    }
    normle =  scone_res$normalized_data[[nom]]
    write.table(x = normle,row.names = TRUE,col.names = TRUE,
                quote = FALSE,sep = "\t",eol = "\n",
                file = paste0(nom_dir,"/norm.txt"))
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
  
  #stop duplicating stdout to log file
  sink(type = "output")
  #sink(type = "message")  -- no, cannot split the message connection
}