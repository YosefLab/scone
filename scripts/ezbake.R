#!/usr/local/bin/Rscript
args <- commandArgs(TRUE)

if(length(args) < 1) {
  args <- c("--help")
}

if("--help" %in% args) {
  cat("
      SCONE Easy-Bake Script: Run Essential SCONE Modules
      
      Arguments:\n
      --expr\tPath to a tab-delimited text file containing all of your gene expression\n
      \t\tscores for all genes and single cells. The file must be a plain text (.txt) file 
      \t\twith with a header row containing the value 'GENE' in the first column, and single 
      \t\tcell names in each successive column.\n
      --qc\tPath to tab-delim txt file containing quality control (QC) matrix (cells in rows, metrics in columns). Must include header.\n
      --bio\t(Optional) Path to txt file containing biological conditions (cell per row).
      \t\tThis condition will be modeled in the Adjustment Step
      \t\tas variation to be preserved. If adjust_bio='no', it will not be used
      \t\tfor normalization, but only for evaluation.\n
      --batch\t(Optional) Path to txt file containing batch conditions (cell per row).
      \t\tThe batch variable will be included in the
      \t\tadjustment model as variation to be removed. If adjust_batch='no', it will
      \t\tnot be used for normalization, but only for evaluation.\n
      --negcon\tPath to txt file containing negative control genes (gene per row). 
      \t\tThese genes are to be used as negative controls for filtering, 
      \t\tnormalization, and evaluation. These genes should be expressed uniformily across
      \t\tthe biological phenomenon of interest.\n
      --out_dir\tOutput directory. Default current wd.\n
      --seed\tRandom seed. Default 112233.\n
      Filtering Arguments:\n
      --filt_genes\tTRUE or FALSE. Should genes be filtered post-sample filtering? Default TRUE.\n
      False-Negative Model Arguments:\n
      --fnr_maxiter\tMaximum number of iterations in EM estimation of expression posteriors. Default 1000.\n
      Normalization Arguments:\n
      --norm_impute\tShould imputation be included in the comparison? 'yes', 'no', or 'force'. Default 'yes'.\n
      --norm_scaling\tScaling options to be included in the Scaling Step, separated by commas. Default 'none,sum,deseq,tmm,uq,fq,detect'.\n
      --norm_rezero\tTRUE or FALSE. Restore imputed zeroes to zero following the Scaling Step? Default TRUE.\n
      --norm_k_max\t(Optional) Maximum number (norm_k_max) of factors of unwanted variation modeled
      \t\tin the Adjustment Step.\n
      --norm_qc_expl\tIn automatic selection of norm_k_max, what
      \t\tfraction of variation must be explained by the first norm_k_max PCs of qc?
      \t\tDefault 0.5. Ignored if norm_k_max is not specified.\n
      --norm_adjust_bio\tIf 'no' it will not be included in the model; if
      \t\t'yes', both models with and without 'bio' will be run; if 'force', only
      \t\tmodels with 'bio' will be run. Default 'yes'.\n
      --norm_adjust_batch\tIf 'no' it will not be modeled in the Adjustment Step;
      \t\tif 'yes', both models with and without 'batch' will be run; if 'force',
      \t\tonly models with 'batch' will be run. Default 'yes'.\n
      Evaluation Arguments:\n
      --eval_dim\t(Optional) The number of principal components to use for evaluation.\n
      --eval_expr_expl\tIn automatic selection of eval_dim, what
      \t\tfraction of variation must be explained by the first eval_dim PCs of expr?
      \t\tDefault 0.1. Ignored if eval_dim is not specified.\n
      --eval_poscon\t(Optional) The genes to be used as positive controls for
      \t\tevaluation. These genes should be expected to change according to the
      \t\tbiological phenomenon of interest.\n
      --eval_max_kclust\tThe max number of clusters (> 1) to be used for pam
      \t\ttightness evaluation. Default 10.\n
      --eval_stratified_pam\tIf TRUE then maximum ASW for PAM_SIL is
      \t\tseparately computed for each biological-cross-batch condition (accepting
      \t\tNAs), and a weighted average is returned as PAM_SIL. Default TRUE.\n
      Additional Arguments:\n
      --report_num\tNumber of top methods to report. Default 13.\n
      --verbose\tVerbosity level (0,1,2): higher level is more verbose. Default 0.\n
      --help\t(Optional) Print this help message.\n
      
      Example:
      ./ezbake.R --help \n\n
      ")
  
  q(save="no")
}

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(is.null(argsL$expr)) {
  stop("expr argument is missing")
}else{
  print("Loading expression matrix...")
  expr_tab = read.table(argsL$expr, sep = "\t",header = TRUE,row.names = NULL)
  expr = as.matrix(expr_tab[,-1])
  genes = as.character(expr_tab$GENE)
  rownames(expr) = genes
  print("Complete")
}

if(is.null(argsL$qc)) {
  stop("qc argument is missing")
}else{
  print("Loading qc matrix...")
  qc = read.table(argsL$qc,sep = "\t",header = TRUE,row.names = NULL)
  rownames(qc) = colnames(expr)
  print("Complete")
}

if(!is.null(argsL$bio)) {
  print("Loading bio data...")
  bio = factor(unlist(read.table(argsL$bio,header = FALSE, sep = "\t",col.names = FALSE)))
  print("Complete")
  names(bio) = colnames(expr)
}else{
  bio = NULL
}

if(!is.null(argsL$batch)) {
  print("Loading batch data...")
  batch = factor(unlist(read.table(argsL$batch,header = FALSE, sep = "\t",col.names = FALSE)))
  print("Complete")
  names(batch) = colnames(expr)
}else{
  batch = NULL
}

if(is.null(argsL$negcon)) {
  stop("negcon argument is missing")
}else{
  print("Loading negative control gene names...")
  negcon = as.character(unlist(read.table(argsL$negcon,header = FALSE, sep = "\t",col.names = FALSE)))
  print("Complete")
}

if(is.null(argsL$out_dir)) {
  out_dir = getwd()
}else{
  out_dir = argsL$out_dir
}

if(is.null(argsL$seed)) {
  seed = 112233
}else{
  seed = as.numeric(argsL$seed)
}

if(is.null(argsL$filt_genes)) {
  filt_genes = TRUE
}else{
  filt_genes = as.logical(argsL$filt_genes)
}

if(is.null(argsL$fnr_maxiter)) {
  fnr_maxiter = 1000
}else{
  fnr_maxiter = as.numeric(argsL$fnr_maxiter)
}

if(is.null(argsL$norm_impute)) {
  norm_impute = "yes"
}else{
  norm_impute = argsL$norm_impute
}

if(is.null(argsL$norm_scaling)) {
  argsL$norm_scaling = "none,sum,deseq,tmm,uq,fq,detect"
}
norm_scaling = strsplit(argsL$norm_scaling,split = ",")[[1]]

if(is.null(argsL$norm_rezero)) {
  norm_rezero = TRUE
}else{
  norm_rezero = as.logical(argsL$norm_rezero)
}

if(is.null(argsL$norm_k_max)) {
  norm_k_max = NULL
}else{
  norm_k_max = as.numeric(argsL$norm_k_max)
}

if(is.null(argsL$norm_qc_expl)) {
  norm_qc_expl = 0.5
}else{
  norm_qc_expl = as.numeric(argsL$norm_qc_expl)
}

if(is.null(argsL$norm_adjust_bio)) {
  norm_adjust_bio = "yes"
}else{
  norm_adjust_bio = argsL$norm_adjust_bio
}

if(is.null(argsL$norm_adjust_batch)) {
  norm_adjust_batch = "yes"
}else{
  norm_adjust_batch = argsL$norm_adjust_batch
}

if(is.null(argsL$eval_dim)) {
  eval_dim = NULL
}else{
  eval_dim = as.numeric(argsL$eval_dim)
}

if(is.null(argsL$eval_expr_expl)) {
  eval_expr_expl = 0.1
}else{
  eval_expr_expl = as.numeric(argsL$eval_expr_expl)
}

if(is.null(argsL$eval_poscon)) {
  eval_poscon = NULL
}else{
  print("Loading positive control gene names...")
  eval_poscon = as.character(unlist(read.table(argsL$eval_poscon,header = FALSE, sep = "\t",col.names = FALSE)))
  print("Complete")
}

if(is.null(argsL$eval_max_kclust)) {
  eval_max_kclust = 10
}else{
  eval_max_kclust = as.numeric(argsL$eval_max_kclust)
}

if(is.null(argsL$eval_stratified_pam)) {
  eval_stratified_pam = TRUE
}else{
  eval_stratified_pam = as.logical(argsL$eval_stratified_pam)
}

if(is.null(argsL$report_num)) {
  report_num = 13
}else{
  report_num = as.numeric(argsL$report_num)
}

if(is.null(argsL$verbose)) {
  verbose = "0"
}else{
  verbose = as.character(argsL$verbose)
}

## ----- SCONE -----
scone::scone_easybake(expr = expr, qc = qc, bio = bio, batch = batch,
               negcon = negcon, verbose = verbose, out_dir = out_dir, seed = seed,
               filt_genes = filt_genes,
               fnr_maxiter = fnr_maxiter,
               norm_impute = norm_impute, norm_scaling = norm_scaling, norm_rezero = norm_rezero,
               norm_k_max = norm_k_max, norm_qc_expl = norm_qc_expl,
               norm_adjust_bio = norm_adjust_bio, norm_adjust_batch = norm_adjust_batch,
               eval_dim = eval_dim, eval_expr_expl = eval_expr_expl, eval_poscon = eval_poscon,
               eval_max_kclust = eval_max_kclust, eval_stratified_pam = eval_stratified_pam,
               report_num = report_num)

q(save="no")
