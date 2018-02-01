pkgname <- "scone"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('scone')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("CLR_FN")
### * CLR_FN

flush(stderr()); flush(stdout())

### Name: CLR_FN
### Title: Centered log-ratio (CLR) normalization wrapper function
### Aliases: CLR_FN

### ** Examples

ei <- matrix(0:20,nrow = 7)
eo <- CLR_FN(ei)




cleanEx()
nameEx("DESEQ_FN")
### * DESEQ_FN

flush(stderr()); flush(stdout())

### Name: DESEQ_FN
### Title: Relative log-expression (RLE; DESeq) scaling normalization
###   wrapper function
### Aliases: DESEQ_FN

### ** Examples

ei <- matrix(0:20,nrow = 7)
eo <- DESEQ_FN(ei)




cleanEx()
nameEx("FQ_FN")
### * FQ_FN

flush(stderr()); flush(stdout())

### Name: FQ_FN
### Title: Full-quantile normalization wrapper function
### Aliases: FQ_FN FQT_FN

### ** Examples

ei <- matrix(0:20,nrow = 7)
eo <- FQ_FN(ei)

ei <- matrix(0:20,nrow = 7)
eo <- FQT_FN(ei)




cleanEx()
nameEx("SCRAN_FN")
### * SCRAN_FN

flush(stderr()); flush(stdout())

### Name: SCRAN_FN
### Title: Simple deconvolution normalization wrapper
### Aliases: SCRAN_FN

### ** Examples

ei <- matrix(0:76,nrow = 7)
eo <- SCRAN_FN(ei)




cleanEx()
nameEx("SUM_FN")
### * SUM_FN

flush(stderr()); flush(stdout())

### Name: SUM_FN
### Title: Sum scaling normalization function
### Aliases: SUM_FN

### ** Examples

ei <- matrix(0:20,nrow = 7)
eo <- SUM_FN(ei)




cleanEx()
nameEx("SconeExperiment-class")
### * SconeExperiment-class

flush(stderr()); flush(stdout())

### Name: SconeExperiment-class
### Title: Class SconeExperiment
### Aliases: SconeExperiment-class SconeExperiment SconeExperiment
###   SconeExperiment,SummarizedExperiment-method
###   SconeExperiment,matrix-method

### ** Examples

set.seed(42)
nrows <- 200
ncols <- 6
counts <- matrix(rpois(nrows * ncols, lambda=10), nrows)
rowdata <- data.frame(poscon=c(rep(TRUE, 10), rep(FALSE, nrows-10)))
coldata <- data.frame(bio=gl(2, 3))
se <- SummarizedExperiment(assays=SimpleList(counts=counts),
                          rowData=rowdata, colData=coldata)

scone1 <- SconeExperiment(assay(se), bio=coldata$bio, poscon=rowdata$poscon)

scone2 <- SconeExperiment(se, which_bio=1L, which_poscon=1L)





cleanEx()
nameEx("TMM_FN")
### * TMM_FN

flush(stderr()); flush(stdout())

### Name: TMM_FN
### Title: Weighted trimmed mean of M-values (TMM) scaling normalization
###   wrapper function
### Aliases: TMM_FN

### ** Examples

ei <- matrix(0:20,nrow = 7)
eo <- TMM_FN(ei)




cleanEx()
nameEx("UQ_FN")
### * UQ_FN

flush(stderr()); flush(stdout())

### Name: UQ_FN
### Title: Upper-quartile (UQ) scaling normalization wrapper function
### Aliases: UQ_FN

### ** Examples

ei <- matrix(0:20,nrow = 7)
eo <- UQ_FN(ei)




cleanEx()
nameEx("biplot_color")
### * biplot_color

flush(stderr()); flush(stdout())

### Name: biplot_color
### Title: Function for biplotting with no point labels and with points
###   color-coded according to a quantitative variable. For example: the
###   rank of normalization performance.
### Aliases: biplot_color

### ** Examples

mat <- matrix(rnorm(1000), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")

pc <- prcomp(mat)

biplot_color(pc, rank(pc$x[,1]))




cleanEx()
nameEx("biplot_interactive")
### * biplot_interactive

flush(stderr()); flush(stdout())

### Name: biplot_interactive
### Title: Interactive biplot
### Aliases: biplot_interactive

### ** Examples

mat <- matrix(rpois(1000, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat)
res <- scone(obj, scaling=list(none=identity,
   uq=UQ_FN, deseq=DESEQ_FN,  fq=FQT_FN),
evaluate=TRUE, k_ruv=0, k_qc=0, eval_kclust=2,
   bpparam = BiocParallel::SerialParam())
## Not run: 
##D biplot_interactive(res)
## End(Not run)




cleanEx()
nameEx("control_genes")
### * control_genes

flush(stderr()); flush(stdout())

### Name: control_genes
### Title: Data: Positive and Negative Control Genes
### Aliases: control_genes cortical_markers housekeeping
###   housekeeping_revised cellcycle_genes

### ** Examples

data(housekeeping)
data(housekeeping_revised)
data(cellcycle_genes)
data(cortical_markers)



cleanEx()
nameEx("estimate_ziber")
### * estimate_ziber

flush(stderr()); flush(stdout())

### Name: estimate_ziber
### Title: Parameter estimation of zero-inflated bernoulli model
### Aliases: estimate_ziber

### ** Examples

mat <- matrix(rpois(1000, lambda = 3), ncol=10)
mat = mat * matrix(1-rbinom(1000, size = 1, prob = .01), ncol=10)
ziber_out = suppressWarnings(estimate_ziber(mat,
   bulk_model = TRUE,
   pos_controls = 1:10))




cleanEx()
nameEx("factor_sample_filter")
### * factor_sample_filter

flush(stderr()); flush(stdout())

### Name: factor_sample_filter
### Title: Factor-based Sample Filtering: Function to filter single-cell
###   RNA-Seq libraries.
### Aliases: factor_sample_filter

### ** Examples

mat <- matrix(rpois(1000, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
qc = as.matrix(cbind(colSums(mat),colSums(mat > 0)))
rownames(qc) = colnames(mat)
colnames(qc) = c("NCOUNTS","NGENES")
mfilt = factor_sample_filter(expr = mat,
    qc, plot = TRUE,qual_select_q_thresh = 1)




cleanEx()
nameEx("fast_estimate_ziber")
### * fast_estimate_ziber

flush(stderr()); flush(stdout())

### Name: fast_estimate_ziber
### Title: Fast parameter estimation of zero-inflated bernoulli model
### Aliases: fast_estimate_ziber

### ** Examples

mat <- matrix(rpois(1000, lambda = 3), ncol=10)
mat = mat * matrix(1-rbinom(1000, size = 1, prob = .01), ncol=10)
ziber_out = suppressWarnings(fast_estimate_ziber(mat,
   bulk_model = TRUE,
   pos_controls = 1:10))




cleanEx()
nameEx("get_bio")
### * get_bio

flush(stderr()); flush(stdout())

### Name: get_bio
### Title: Get Factor of Biological Conditions and Batch
### Aliases: get_bio get_batch get_bio,SconeExperiment-method
###   get_batch,SconeExperiment-method get_batch
###   get_bio,SconeExperiment-method get_batch,SconeExperiment-method

### ** Examples

set.seed(42)
mat <- matrix(rpois(500, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat, bio = factor(rep(c(1,2),each = 5)),
           batch = factor(rep(c(1,2),times = 5)))
bio = get_bio(obj)
batch = get_batch(obj)




cleanEx()
nameEx("get_design")
### * get_design

flush(stderr()); flush(stdout())

### Name: get_design
### Title: Retrieve Design Matrix
### Aliases: get_design get_design,SconeExperiment,character-method
###   get_design,SconeExperiment,numeric-method

### ** Examples

set.seed(42)
mat <- matrix(rpois(500, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat, bio = factor(rep(c(1,2),each = 5)),
           batch = factor(rep(c(1,2),times = 5)))
res <- scone(obj, scaling=list(none=identity, uq=UQ_FN),
           evaluate=TRUE, k_ruv=0, k_qc=0, 
           adjust_batch = "yes", adjust_bio = "yes",
           eval_kclust=2, bpparam = BiocParallel::SerialParam())
design_top = get_design(res,1)




cleanEx()
nameEx("get_negconruv")
### * get_negconruv

flush(stderr()); flush(stdout())

### Name: get_negconruv
### Title: Get Negative and Positive Controls
### Aliases: get_negconruv get_negconeval get_poscon
###   get_negconruv,SconeExperiment-method
###   get_negconeval,SconeExperiment-method
###   get_poscon,SconeExperiment-method get_negconeval get_poscon
###   get_negconruv,SconeExperiment-method
###   get_negconeval,SconeExperiment-method
###   get_poscon,SconeExperiment-method

### ** Examples

set.seed(42)
mat <- matrix(rpois(500, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat,negcon_ruv = 1:50 %in% 1:10,
           negcon_eval = 1:50 %in% 11:20,
           poscon = 1:50 %in% 21:30)
negcon_ruv = get_negconruv(obj)
negcon_eval = get_negconeval(obj)
poscon = get_poscon(obj)




cleanEx()
nameEx("get_normalized")
### * get_normalized

flush(stderr()); flush(stdout())

### Name: get_normalized
### Title: Retrieve Normalized Matrix
### Aliases: get_normalized get_normalized,SconeExperiment,character-method
###   get_normalized,SconeExperiment,numeric-method

### ** Examples

set.seed(42)
mat <- matrix(rpois(500, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat)
res <- scone(obj, scaling=list(none=identity, uq=UQ_FN),
           evaluate=TRUE, k_ruv=0, k_qc=0, 
           eval_kclust=2, bpparam = BiocParallel::SerialParam())
top_norm = get_normalized(res,1)
           




cleanEx()
nameEx("get_params")
### * get_params

flush(stderr()); flush(stdout())

### Name: get_params
### Title: Extract scone parameters
### Aliases: get_params get_params,SconeExperiment-method
###   get_params,SconeExperiment-method

### ** Examples

set.seed(42)
mat <- matrix(rpois(500, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat)
res <- scone(obj, scaling=list(none=identity, uq=UQ_FN),
           run = FALSE, k_ruv=0, k_qc=0, eval_kclust=2)
params = get_params(res)




cleanEx()
nameEx("get_qc")
### * get_qc

flush(stderr()); flush(stdout())

### Name: get_qc
### Title: Get Quality Control Matrix
### Aliases: get_qc get_qc,SconeExperiment-method
###   get_qc,SconeExperiment-method

### ** Examples

set.seed(42)
mat <- matrix(rpois(500, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat,
         qc = cbind(colSums(mat),colSums(mat > 0)))
qc = get_qc(obj)




cleanEx()
nameEx("get_scores")
### * get_scores

flush(stderr()); flush(stdout())

### Name: get_scores
### Title: Extract scone scores
### Aliases: get_scores get_scores,SconeExperiment-method get_score_ranks
###   get_score_ranks,SconeExperiment-method get_score_ranks
###   get_scores,SconeExperiment-method
###   get_score_ranks,SconeExperiment-method

### ** Examples

set.seed(42)
mat <- matrix(rpois(500, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat)
res <- scone(obj, scaling=list(none=identity, uq=UQ_FN),
           evaluate=TRUE, k_ruv=0, k_qc=0, 
           eval_kclust=2, bpparam = BiocParallel::SerialParam())
scores = get_scores(res)
score_ranks = get_score_ranks(res)




cleanEx()
nameEx("impute_expectation")
### * impute_expectation

flush(stderr()); flush(stdout())

### Name: impute_expectation
### Title: Imputation of zero abundance based on general zero-inflated
###   model
### Aliases: impute_expectation

### ** Examples

mat <- matrix(rpois(1000, lambda = 3), ncol=10)
mat = mat * matrix(1-rbinom(1000, size = 1, prob = .01), ncol=10)

mu = matrix(rep(3/ppois(0,lambda = 3,lower.tail = FALSE),1000),ncol = 10)

p_false = 1 / ( 1 + ppois(0, lambda = 3, lower.tail = TRUE ) / 
    (0.01 * ppois(0, lambda = 3, lower.tail = FALSE) ) )

p_nodrop = matrix(rep(1-p_false,1000),ncol = 10)
p_nodrop[mat > 0] = 1

impute_args = list()
impute_args = list(mu = mu, p_nodrop = p_nodrop)

imat = impute_expectation(mat,impute_args = impute_args)




cleanEx()
nameEx("impute_null")
### * impute_null

flush(stderr()); flush(stdout())

### Name: impute_null
### Title: Null or no-op imputation
### Aliases: impute_null

### ** Examples

mat <- matrix(rpois(1000, lambda = 5), ncol=10)
imat = impute_null(mat)




cleanEx()
nameEx("lm_adjust")
### * lm_adjust

flush(stderr()); flush(stdout())

### Name: lm_adjust
### Title: Linear Adjustment Normalization
### Aliases: lm_adjust

### ** Examples


set.seed(141)
bio = as.factor(rep(c(1,2),each = 2))
batch = as.factor(rep(c(1,2),2))
design_mat = make_design(bio,batch, W = NULL)

log_expr = matrix(rnorm(20),ncol = 4)
adjusted_log_expr = lm_adjust(log_expr = log_expr,
  design_mat = design_mat,
  batch = batch)




cleanEx()
nameEx("make_design")
### * make_design

flush(stderr()); flush(stdout())

### Name: make_design
### Title: Make a Design Matrix
### Aliases: make_design

### ** Examples


bio = as.factor(rep(c(1,2),each = 2))
batch = as.factor(rep(c(1,2),2))
design_mat = make_design(bio,batch, W = NULL)




cleanEx()
nameEx("metric_sample_filter")
### * metric_sample_filter

flush(stderr()); flush(stdout())

### Name: metric_sample_filter
### Title: Metric-based Sample Filtering: Function to filter single-cell
###   RNA-Seq libraries.
### Aliases: metric_sample_filter

### ** Examples

mat <- matrix(rpois(1000, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
qc = as.matrix(cbind(colSums(mat),colSums(mat > 0)))
rownames(qc) = colnames(mat)
colnames(qc) = c("NCOUNTS","NGENES")
mfilt = metric_sample_filter(expr = mat,nreads = qc[,"NCOUNTS"],
   plot = TRUE, hard_nreads = 0)




cleanEx()
nameEx("scone")
### * scone

flush(stderr()); flush(stdout())

### Name: scone
### Title: Normalize Expression Data and Evaluate Normalization Performance
### Aliases: scone scone scone,SconeExperiment-method

### ** Examples


mat <- matrix(rpois(1000, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat)
no_results <- scone(obj, scaling=list(none=identity,
           uq=UQ_FN, deseq=DESEQ_FN),
           run=FALSE, k_ruv=0, k_qc=0, eval_kclust=2)
           
results <- scone(obj, scaling=list(none=identity,
           uq=UQ_FN, deseq=DESEQ_FN),
           run=TRUE, k_ruv=0, k_qc=0, eval_kclust=2,
           bpparam = BiocParallel::SerialParam())
           
results_in_memory <- scone(obj, scaling=list(none=identity,
           uq=UQ_FN, deseq=DESEQ_FN),
           k_ruv=0, k_qc=0, eval_kclust=2,
           return_norm = "in_memory",
           bpparam = BiocParallel::SerialParam())




cleanEx()
nameEx("sconeReport")
### * sconeReport

flush(stderr()); flush(stdout())

### Name: sconeReport
### Title: SCONE Report Browser: Browse Evaluation of Normalization
###   Performance
### Aliases: sconeReport

### ** Examples

set.seed(101)
mat <- matrix(rpois(1000, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat)
res <- scone(obj, scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
           evaluate=TRUE, k_ruv=0, k_qc=0, eval_kclust=2,
           bpparam = BiocParallel::SerialParam())
qc = as.matrix(cbind(colSums(mat),colSums(mat > 0)))
rownames(qc) = colnames(mat)
colnames(qc) = c("NCOUNTS","NGENES")
## Not run: 
##D sconeReport(res,rownames(get_params(res)), qc = qc)
## End(Not run)




cleanEx()
nameEx("scone_easybake")
### * scone_easybake

flush(stderr()); flush(stdout())

### Name: scone_easybake
### Title: Wrapper for Running Essential SCONE Modules
### Aliases: scone_easybake

### ** Examples

set.seed(101)
mat <- matrix(rpois(1000, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat)
res <- scone(obj, scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
           evaluate=TRUE, k_ruv=0, k_qc=0, eval_kclust=2, 
           bpparam = BiocParallel::SerialParam())
qc = as.matrix(cbind(colSums(mat),colSums(mat > 0)))
rownames(qc) = colnames(mat)
colnames(qc) = c("NREADS","RALIGN")
## Not run: 
##D scone_easybake(mat, qc = as.data.frame(qc), verbose = "2", 
##D    norm_adjust_bio= "no",
##D    norm_adjust_batch= "no", norm_k_max = 0,
##D    fnr_maxiter = 0, filt_cells=FALSE, filt_genes=FALSE,
##D    eval_stratified_pam = FALSE,
##D    out_dir="~/scone_out")
## End(Not run)



cleanEx()
nameEx("score_matrix")
### * score_matrix

flush(stderr()); flush(stdout())

### Name: score_matrix
### Title: SCONE Evaluation: Evaluate an Expression Matrix
### Aliases: score_matrix

### ** Examples


set.seed(141)
bio = as.factor(rep(c(1,2),each = 2))
batch = as.factor(rep(c(1,2),2))
log_expr = matrix(rnorm(20),ncol = 4)

scone_metrics = score_matrix(log_expr,
   bio = bio, batch = batch,
   eval_kclust = 2, is_log = TRUE)




cleanEx()
nameEx("select_methods")
### * select_methods

flush(stderr()); flush(stdout())

### Name: select_methods
### Title: Get a subset of normalizations from a SconeExperiment object
### Aliases: select_methods select_methods,SconeExperiment,character-method
###   select_methods,SconeExperiment,numeric-method

### ** Examples

set.seed(42)
mat <- matrix(rpois(500, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat)
res <- scone(obj, scaling=list(none=identity, uq=UQ_FN),
           evaluate=TRUE, k_ruv=0, k_qc=0, 
           eval_kclust=2, bpparam = BiocParallel::SerialParam())
select_res = select_methods(res,1:2)




cleanEx()
nameEx("simple_FNR_params")
### * simple_FNR_params

flush(stderr()); flush(stdout())

### Name: simple_FNR_params
### Title: Fit Simple False-Negative Model
### Aliases: simple_FNR_params

### ** Examples

mat <- matrix(rpois(1000, lambda = 3), ncol=10)
mat = mat * matrix(1-rbinom(1000, size = 1, prob = .01), ncol=10)
fnr_out = simple_FNR_params(mat,pos_controls = 1:10)




cleanEx()
nameEx("subsample_cells")
### * subsample_cells

flush(stderr()); flush(stdout())

### Name: subsample_cells
### Title: Function to subsample Scone object by subsampling cells
### Aliases: subsample_cells

### ** Examples

# Create a scone object
mat <- matrix(rpois(1000, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat)
res <- scone(obj, scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
           evaluate=TRUE, k_ruv=0, k_qc=0, eval_kclust=2,
           bpparam = BiocParallel::SerialParam())
qc = as.matrix(cbind(colSums(mat),colSums(mat > 0)))
rownames(qc) = colnames(mat)
colnames(qc) = c("NCOUNTS","NGENES")

# Subsample the scone object
subsampled_res <- subsample_cells(res,percent=50, at_bio = FALSE, seed = 100, verbose = TRUE)





cleanEx()
nameEx("subsample_cells_with_min_bio")
### * subsample_cells_with_min_bio

flush(stderr()); flush(stdout())

### Name: subsample_cells_with_min_bio
### Title: Internal function to subsample Scone object by subsampling cells
###   by size of smallest bio group
### Aliases: subsample_cells_with_min_bio

### ** Examples

# Create a scone object
mat <- matrix(rpois(1000, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat)
res <- scone(obj, scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
           evaluate=TRUE, k_ruv=0, k_qc=0, eval_kclust=2,
           bpparam = BiocParallel::SerialParam())
qc = as.matrix(cbind(colSums(mat),colSums(mat > 0)))
rownames(qc) = colnames(mat)
colnames(qc) = c("NCOUNTS","NGENES")
res$bio <- c('a', 'a', 'a', 'a', 'a', 'a', 'a', 'b', 'b', 'b')

# Subsample the scone object
subsampled_res <- subsample_cells_with_min_bio(res, seed = 100, verbose = TRUE)





cleanEx()
nameEx("subsample_genes")
### * subsample_genes

flush(stderr()); flush(stdout())

### Name: subsample_genes
### Title: Function to subsample Scone object by subsampling genes
### Aliases: subsample_genes

### ** Examples

# Create a scone object
mat <- matrix(rpois(1000, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat)
res <- scone(obj, scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
           evaluate=TRUE, k_ruv=0, k_qc=0, eval_kclust=2,
           bpparam = BiocParallel::SerialParam())
qc = as.matrix(cbind(colSums(mat),colSums(mat > 0)))
rownames(qc) = colnames(mat)
colnames(qc) = c("NCOUNTS","NGENES")

# Subsample the scone object
subsampled_res <- subsample_genes(res, percent=50, keep_all_control = TRUE, seed = 100, verbose = TRUE)





cleanEx()
nameEx("subsample_scone")
### * subsample_scone

flush(stderr()); flush(stdout())

### Name: subsample_scone
### Title: Wrapper Function to subsample Scone object
### Aliases: subsample_scone

### ** Examples

# Create a scone object
mat <- matrix(rpois(1000, lambda = 5), ncol=10)
colnames(mat) <- paste("X", 1:ncol(mat), sep="")
obj <- SconeExperiment(mat)
scone_object <- scone(obj, scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
           evaluate=TRUE, k_ruv=0, k_qc=0, eval_kclust=2,
           bpparam = BiocParallel::SerialParam())
qc = as.matrix(cbind(colSums(mat),colSums(mat > 0)))
rownames(qc) = colnames(mat)
colnames(qc) = c("NCOUNTS","NGENES")

scone_object$bio <- c('a', 'a', 'a', 'a', 'a', 'a', 'a', 'b', 'b', 'b')

# Subsample the scone object
subsampled_res <- subsample_scone(scone_object, subsample_gene_level = 50, subsample_cell_level = 100, 
                                     cells_first=TRUE, at_bio = TRUE,
                                     keep_all_control=TRUE, seed = 100, verbose = TRUE)





### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
