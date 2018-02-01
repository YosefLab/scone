## ----options, include=FALSE, cache=FALSE, results='hide', message=FALSE----

knitr::opts_chunk$set(fig.align="center", cache=FALSE,error=FALSE,
                      fig.width=6,fig.height=6,autodep=TRUE,
                      out.width="600px", out.height="600px",
                      results="markup", echo=TRUE, eval=TRUE)

options(getClass.msg=FALSE)

set.seed(6473) ## for reproducibility

library(scone)
library(RColorBrewer)


## ----datain--------------------------------------------------------------

library(scRNAseq)

## ----- Load Example Data -----
data(fluidigm)

# Set assay to RSEM estimated counts
assay(fluidigm) = assays(fluidigm)$rsem_counts


## ----showqc--------------------------------------------------------------

## ----- List all QC fields -----

# List all qc fields (accessible via colData())
metadata(fluidigm)$which_qc


## ----biocoverage---------------------------------------------------------

# Joint distribution of "biological condition"" and "coverage type""
table(colData(fluidigm)$Coverage_Type,
      colData(fluidigm)$Biological_Condition)


## ----prefilter-----------------------------------------------------------

# Preliminary Sample Filtering: High-Coverage Only
is_select = colData(fluidigm)$Coverage_Type == "High"
fluidigm = fluidigm[,is_select]

# Retain only detected transcripts
fluidigm = fluidigm[which(apply(assay(fluidigm) > 0,1,any)),]


## ----ralign--------------------------------------------------------------

# Define a color scheme
cc <- c(brewer.pal(9, "Set1"))

# One batch per Biological Condition
batch = factor(colData(fluidigm)$Biological_Condition)

# Alignment Quality Metrics
qc = colData(fluidigm)[,metadata(fluidigm)$which_qc]

# Barplot of read proportion mapping to human transcriptome
ralign = qc$RALIGN
o = order(ralign)[order(batch[order(ralign)])] # Order by batch, then value

barplot(ralign[o], col=cc[batch][o], 
        border=cc[batch][o], main="Percentage of reads mapped")
legend("bottomleft", legend=levels(batch), fill=cc,cex=0.4)


## ----nreads--------------------------------------------------------------

# Barplot of total read number
nreads = qc$NREADS
o = order(nreads)[order(batch[order(nreads)])] # Order by batch, then value

barplot(nreads[o], col=cc[batch][o], 
        border=cc[batch][o], main="Total number of reads")
legend("topright", legend=levels(batch), fill=cc, cex=0.4)


## ----qpc-----------------------------------------------------------------

## ----- PCA of QC matrix -----
qpc = prcomp(qc,center = TRUE,scale. = TRUE)
barplot((qpc$sdev^2)/sum(qpc$sdev^2), border="gray", 
        xlab="PC", ylab="Proportion of Variance", main="Quality PCA")


## ----qpc_view------------------------------------------------------------

# Barplot of PC1 of the QC matrix
qc1 = as.vector(qpc$x[,1])
o = order(qc1)[order(batch[order(qc1)])]

barplot(qc1[o], col=cc[batch][o], 
        border=cc[batch][o], main="Quality PC1")
legend("bottomright", legend=levels(batch), 
       fill=cc, cex=0.8)


## ----fnr_fit-------------------------------------------------------------

# Extract Housekeeping Genes
data(housekeeping)
hk = intersect(housekeeping$V1,rownames(assay(fluidigm)))

# Mean log10(x+1) expression
mu_obs = rowMeans(log10(assay(fluidigm)[hk,]+1))

# Assumed False Negatives
drop_outs = assay(fluidigm)[hk,] == 0

# Logistic Regression Model of Failure
ref.glms = list()
for (si in 1:dim(drop_outs)[2]){
  fit = glm(cbind(drop_outs[,si],1 - drop_outs[,si]) ~ mu_obs,
            family=binomial(logit))
  ref.glms[[si]] = fit$coefficients
}


## ----fnr_vis,fig.width=8,fig.height=4,out.width="800px",out.height="400px"----

par(mfrow=c(1,2))

# Plot Failure Curves and Calculate AUC
plot(NULL, main = "False Negative Rate Curves",
     ylim = c(0,1),xlim = c(0,6), 
     ylab = "Failure Probability", xlab = "Mean log10 Expression")
x = (0:60)/10
AUC = NULL
for(si in 1:ncol(assay(fluidigm))){
  y = 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1)
  AUC[si] = sum(y)/10
  lines(x, 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1),
        type = 'l', lwd = 2, col = cc[batch][si])
}

# Barplot of FNR AUC
o = order(AUC)[order(batch[order(AUC)])]

barplot(AUC[o], col=cc[batch][o], border=cc[batch][o], main="FNR AUC")
legend("topright", legend=levels(batch), fill=cc, cex=0.4)


## ----metric_sample_filter, fig.height= 10,out.height="1000px"------------

# Initial Gene Filtering: 
# Select "common" transcripts based on proportional criteria.
num_reads = quantile(assay(fluidigm)[assay(fluidigm) > 0])[4]
num_cells = 0.25*ncol(fluidigm)
is_common = rowSums(assay(fluidigm) >= num_reads ) >= num_cells

# Metric-based Filtering
mfilt = metric_sample_filter(assay(fluidigm),
                             nreads = colData(fluidigm)$NREADS,
                             ralign = colData(fluidigm)$RALIGN,
                             gene_filter = is_common,
                             pos_controls = rownames(fluidigm) %in% hk,

                             zcut = 3, mixture = FALSE,
                             plot = TRUE)

# Simplify to a single logical
mfilt = !apply(simplify2array(mfilt[!is.na(mfilt)]),1,any)


## ----thresh,fig.width= 6, fig.height= 4, out.width="600px",out.height="400px"----

hist(qc$RALIGN, breaks = 0:100)
# Hard threshold
abline(v = 15, col = "yellow", lwd = 2)
# 3 (zcut) standard deviations below the mean ralign value
abline(v = mean(qc$RALIGN) - 3*sd(qc$RALIGN), col = "green", lwd = 2)
# 3 (zcut) MADs below the median ralign value
abline(v = median(qc$RALIGN) - 3*mad(qc$RALIGN), col = "red", lwd = 2)
# Sufficient threshold
abline(v = NULL, col = "grey", lwd = 2)

# Final threshold is the minimum of 
# 1) the sufficient threshold and 
# 2) the max of all others
thresh = min(NULL,
             max(c(15,mean(qc$RALIGN) - 3*sd(qc$RALIGN),
                   median(qc$RALIGN) - 3*mad(qc$RALIGN))))
abline(v = thresh, col = "blue", lwd = 2, lty = 2)

legend("topleft",legend = c("Hard","SD","MAD","Sufficient","Final"),
       lwd = 2, col = c("yellow","green","red","grey","blue"),
       lty = c(1,1,1,1,2), cex = .5)


## ----filterCount---------------------------------------------------------

goodDat = fluidigm[,mfilt]

# Final Gene Filtering: Highly expressed in at least 5 cells
num_reads = quantile(assay(fluidigm)[assay(fluidigm) > 0])[4]
num_cells = 5
is_quality = rowSums(assay(fluidigm) >= num_reads ) >= num_cells



## ----scone_init----------------------------------------------------------


# Expression Data (Required)
expr = assay(goodDat)[is_quality,]

# Biological Origin - Variation to be preserved (Optional)
bio = factor(colData(goodDat)$Biological_Condition)

# Processed Alignment Metrics - Variation to be removed (Optional)
qc = colData(goodDat)[,metadata(goodDat)$which_qc]
ppq = scale(qc[,apply(qc,2,sd) > 0],center = TRUE,scale = TRUE)

# Positive Control Genes - Prior knowledge of DE (Optional)
poscon = intersect(rownames(expr),strsplit(paste0("ALS2, CDK5R1, CYFIP1,",
                                                  " DPYSL5, FEZ1, FEZ2, ",
                                                  "MAPT, MDGA1, NRCAM, ",
                                                  "NRP1, NRXN1, OPHN1, ",
                                                  "OTX2, PARD6B, PPT1, ",
                                                  "ROBO1, ROBO2, RTN1, ",
                                                  "RTN4, SEMA4F, SIAH1, ",
                                                  "SLIT2, SMARCA1, THY1, ",
                                                  "TRAPPC4, UBB, YWHAG, ",
                                                  "YWHAH"),split = ", ")[[1]])

# Negative Control Genes - Uniformly expressed transcripts (Optional)
negcon = intersect(rownames(expr),hk)

# Creating a SconeExperiment Object
my_scone <- SconeExperiment(expr,
                qc=ppq, bio = bio,
                negcon_ruv = rownames(expr) %in% negcon,
                poscon = rownames(expr) %in% poscon
)


## ----scone_in2-----------------------------------------------------------

## ----- User-defined function: Dividing by number of detected genes -----

EFF_FN = function (ei)
{
  sums = colSums(ei > 0)
  eo = t(t(ei)*sums/mean(sums))
  return(eo)
}

## ----- Scaling Argument -----

scaling=list(none=identity, # Identity - do nothing

             eff = EFF_FN, # User-defined function

             sum = SUM_FN, # SCONE library wrappers...
             tmm = TMM_FN, 
             uq = UQ_FN,
             fq = FQT_FN,
             deseq = DESEQ_FN)


## ----scone_in3, eval=FALSE-----------------------------------------------
#  
#  # Simple FNR model estimation with SCONE::estimate_ziber
#  fnr_out = estimate_ziber(x = expr, bulk_model = TRUE,
#                           pos_controls = rownames(expr) %in% hk,
#                           maxiter = 10000)
#  
#  ## ----- Imputation List Argument -----
#  imputation=list(none=impute_null, # No imputation
#                  expect=impute_expectation) # Replace zeroes
#  
#  ## ----- Imputation Function Arguments -----
#  # accessible by functions in imputation list argument
#  impute_args = list(p_nodrop = fnr_out$p_nodrop, mu = exp(fnr_out$Alpha[1,]))
#  
#  my_scone <- scone(my_scone,
#                  imputation = imputation, impute_args = impute_args,
#                  scaling=scaling,
#                  k_qc=3, k_ruv = 3,
#                  adjust_bio="no",
#                  run=FALSE)

## ----scone_params--------------------------------------------------------

my_scone <- scone(my_scone,
                scaling=scaling,
                k_qc=3, k_ruv = 3,
                adjust_bio="no",
                run=FALSE)

head(get_params(my_scone))


## ----scone_params_view---------------------------------------------------

apply(get_params(my_scone),2,unique)


## ----scone_params_filt, eval=FALSE---------------------------------------
#  
#  is_screened = ((get_params(my_scone)$imputation_method == "expect") &
#                   (get_params(my_scone)$scaling_method %in% c("none",
#                                                               "eff")))
#  
#  my_scone = select_methods(my_scone,
#                            rownames(get_params(my_scone))[!is_screened ])
#  

## ----scone_run-----------------------------------------------------------

BiocParallel::register(
  BiocParallel::SerialParam()
) # Register BiocParallel Serial Execution

my_scone <- scone(my_scone,
                  scaling=scaling,
                  run=TRUE,
                  eval_kclust = 2:6,stratified_pam = TRUE,
                  return_norm = "in_memory",
                  zero = "postadjust")


## ----scone_view1---------------------------------------------------------

# View Metric Scores
head(get_scores(my_scone))

# View Mean Score Rank
head(get_score_ranks(my_scone))

# Extract normalized data from top method
out_norm = get_normalized(my_scone,
                          method = rownames(get_params(my_scone))[1])


## ----biplot_color--------------------------------------------------------

pc_obj = prcomp(apply(t(get_scores(my_scone)),1,rank),
                center = TRUE,scale = FALSE)
bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)


## ----biplot_color4-------------------------------------------------------

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


## ----sconeReport, eval=FALSE---------------------------------------------
#  
#  # Methods to consider
#  scone_methods = c(rownames(get_params(my_scone))[1:12],
#                    "none,none,no_uv,no_bio,no_batch")
#  
#  # Shiny app
#  sconeReport(my_scone,methods = scone_methods,
#              qc = ppq,
#              bio = bio,
#              negcon = negcon, poscon = poscon)
#  

## ----session-------------------------------------------------------------
sessionInfo()

