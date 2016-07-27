#' Wrapper for Running Essential SCONE Modules
#'
#' @param expr matrix. The expression data matrix (genes in rows, cells in columns).
#' @param qc matrix. The quality control (QC) matrix (cells in rows, metrics in columns)
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
#' @param verbose character. Verbosity level.
#'   Default "0".
#' @param out_dir character. Output directory. 
#'   Default getwd().
#' @param seed numeric. Random seed. 
#'   Default 112233.
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
#'   Default c("sum", "deseq", "tmm", "uq", "fq"). See details.
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
#'   Default 0.2. Ignored if eval_dim is not NULL.
#' @param eval_poscon character. The genes to be used as positive controls for
#'   evaluation. These genes should be expected to change according to the
#'   biological phenomenon of interest.
#' @param eval_kclust numeric. The number of clusters (> 1) to be used for pam
#'   tightness evaluation. If an array of integers, largest average silhouette
#'   width (tightness) will be reported. If NULL, tightness will be returned NA.
#' @param eval_stratified_pam logical. If TRUE then maximum ASW for PAM_SIL is
#'   separately computed for each biological-cross-batch condition (accepting
#'   NAs), and a weighted average is returned as PAM_SIL.
#'   Default TRUE.
#'
#' @export
#' @importFrom RColorBrewer brewer.pal
#' @importFrom boot inv.logit
#' @importFrom hexbin plot hexbin
#'
#' @details "ADD DESCRIPTION"
#'
#' @return Directory structure "ADD DESCRIPTION"
#'

scone_wrapper <- function(expr, qc,
                  bio=NULL, batch=NULL, negcon=NULL,
                  verbose=c("0","1","2"), out_dir=getwd(), seed = 112233,
                  
                  filt_genes=TRUE,
                  
                  fnr_maxiter=1000,
                  
                  norm_impute=c("yes", "no", "force"),
                  norm_scaling=c("sum", "deseq", "tmm", "uq", "fq"),
                  norm_rezero=TRUE, norm_k_max=NULL, norm_qc_expl=0.5,
                  norm_adjust_bio=c("yes", "no", "force"),
                  norm_adjust_batch=c("yes", "no", "force"), 
                  
                  eval_dim=NULL, eval_expr_expl=0.2,
                  eval_poscon=NULL, eval_kclust = 2:10, eval_stratified_pam = TRUE) {

  verbose = match.arg(verbose)
  verbose = as.numeric(verbose)
  
  printf <- function(...) cat(sprintf(...))
  set.seed(seed)
  
  if(!is.matrix(expr)){stop("expr must be a matrix")}
  
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
  
  # Initial Gene Filtering: Select "common" transcripts.
  thresh_fail = quantile(expr[expr > 0], 0.2) ###Make argument###
  num_fail = 10 ###Make argument###
  init.gf.vec = rowSums(expr > thresh_fail) > num_fail
  if(verbose > 0){printf("Data Filtering Module: Initial filter - Kept only %d genes expressed in more than %.2f units in more than %d cells, excluded %d genes\n", sum(init.gf.vec), thresh_fail, num_fail, sum(!init.gf.vec))}
  
  # Initial-Filtering Negative Controls
  small_negcon = intersect(negcon,rownames(expr)[init.gf.vec])
  if(verbose > 0){printf("Data Filtering Module: %d negative control genes to be used for sample filtering\n", length(small_negcon))}
  
  # Metric-based Filtering and Report
  pdf(paste0(misc_dir,"/filter_report.pdf"))
  mfilt = metric_sample_filter(expr, mixture = TRUE, plot = TRUE, hist_breaks = 100,
                               zcut = 4, ###Make argument###
                               pos_controls = rownames(ct) %in% small_negcon,
                               gene_filter = init.gf.vec,
                               nreads = qc$NREADS,ralign = qc$RALIGN) ###Make argument###
  dev.off()
  mfilt = !apply(simplify2array(mfilt[!is.na(mfilt)]),1,any)
  if(verbose > 0){printf("Data Filtering Module: %d samples to be used in downstream analysis, excluded %d samples\n", sum(mfilt),  sum(!mfilt))}
  
  # Implement sample filter
  expr = expr[,mfilt]
  qc = qc[mfilt,]
  batch = droplevels(batch[mfilt])
  bio = droplevels(bio[mfilt])

  if(filt_genes){
    thresh_fail = quantile(expr[expr > 0], 0.2) ###Make argument###
    num_fail = 10 ###Make argument###
    final.gf.vec = rowSums(expr > thresh_fail) > num_fail
    if(verbose > 0){printf("Data Filtering Module: Kept only %d genes expressed in more than %.2f units in more than %d cells, excluded %d genes\n", sum(final.gf.vec), thresh_fail, num_fail, sum(!final.gf.vec))}
  
    # Implement gene filter
    expr = expr[final.gf.vec,]
    if(!is.null(negcon)){
      negcon = negcon[negcon %in% rownames(expr)]
      if(verbose > 0){printf("Data Filtering Module: Kept only %d negative control genes\n", length(negcon))}
    }
    if(!is.null(eval_poscon)){
      eval_poscon = eval_poscon[eval_poscon %in% rownames(expr)]
      if(verbose > 0){printf("Data Filtering Module: Kept only %d positive control genes\n", length(eval_poscon))}
    }
  }
  
  if(verbose > 0){printf("Data Filtering Module: COMPLETE\n\n")}
  ## ------ False-Negative Rate Inference Module ------
  
  if(verbose > 0){printf("False-Negative Rate Inference Module: Estimating posterior expression probability...\n")}
  fnr_out = estimate_ziber(x = expr, 
                           bulk_model = TRUE, ###Make argument###
                           pos_controls = rownames(expr) %in% hk,
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
  
  ui <- fluidPage(
    h1("Example app"),
    sidebarLayout(
      sidebarPanel(
        actionButton("qtBtn", "Quit")
      ),
      mainPanel(
      )
    )
  )

  server <- function(input, output, session) {
    
    observeEvent(input$qtBtn,{
      if(input$qtBtn > 0){
        stopApp("hello world")
      }
    })
  }
  
  print(runApp(list(ui = ui, server = server)))

  sessionInfo()
  
}