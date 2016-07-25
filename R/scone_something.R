library(Biobase) #for AnnotatedDataFrame if scone is not loaded
library(scone)
library(RColorBrewer)
library(boot) #for logit

#A wrapper for filtering and running scone
#currently supports only one-variable bio and batch
scone_something = function(selected_eSet, housekeeping_list,
                           batch, bio = NULL, de = NULL,
                           adjust_bio = "no", adjust_batch = "yes", rezero = TRUE,
                           out_dir = getwd(), stratified_pam = TRUE,
                           selected_expression_in_eSet = "counts_table",
                           results_image_filename = "sconeResults.Rda",
                           plots_filename = "scone_plots.pdf") {

printf <- function(...) cat(sprintf(...))
set.seed(112233) ## reproducibility
if (!file.exists(out_dir)) {
  dir.create(out_dir, recursive=TRUE)
}
pdf(file.path(out_dir, plots_filename), onefile=TRUE)

#keep only protein coding genes
keepGenes = featureData(selected_eSet)$Transcript_Type == "protein_coding"
selected_eSet = selected_eSet[keepGenes,]
printf("Removed %d non-protein coding genes, %d protein-coding genes left\n", sum(!keepGenes), sum(keepGenes))


#don't use exprs because we have several units (e.g., counts, tpm, fpkm) for each eSet
#instead, dynamically access the assayData with the user's argument
selected_unit_exprs = get(selected_expression_in_eSet, envir = assayData(selected_eSet))
selected_unit_exprs[is.na(selected_unit_exprs)] = 0
E = aggregate(selected_unit_exprs, by=list(featureData(selected_eSet)$Gene_Symbol),
              FUN=max, na.rm=TRUE)
rownames(E) = E[,1]
E = E[, -1]
E = data.matrix(E)
E[E == -Inf] = 0 #cases in which the max had no non-NA option - should not happen because we turned NAs into 0's already
stopifnot(all(is.finite(E)))
printf("collapsed %d transcripts into %d genes by gene symbols\n", nrow(selected_unit_exprs), nrow(E))

#now to remove duplicated features from the featuresData (we know that the same gene
#symbol will have the same type, so no harm in just removing duplicate rows) and then
#reordering it according to the new expression matrix
featureDat = featureData(selected_eSet)
featureDat = featureDat[!duplicated(featureDat$Gene_Symbol)]
featureDat = featureDat[rownames(E), ]

# Expression set object to make filtering easier
chosen_units_eSet = new("ExpressionSet", exprs = E, featureData = featureDat, 
                        protocolData = protocolData(selected_eSet), phenoData = phenoData(selected_eSet))

# Remove failed samples:
is.failed = apply(is.na(exprs(chosen_units_eSet)) | (exprs(chosen_units_eSet) == 0), 2, all) |
  apply(is.na(pData(protocolData(chosen_units_eSet))), 1, all)
chosen_units_eSet = chosen_units_eSet[,!is.failed ]
printf("Removed %d failed samples\n", sum(is.failed))


## ----- Pre-Filtering of Transcripts: ----
# Select only coding transcripts and TCR segments
is.expressed.sc = rowMeans(exprs(chosen_units_eSet)) > 0
chosen_units_eSet = chosen_units_eSet[is.expressed.sc, ]
printf("Removed %d undetected genes, %d left\n",sum(!is.expressed.sc), sum(is.expressed.sc))


# Params for gene filtering
# turns out that the filtering of the Th17 dataset was pretty robust to the choice of the quantile
thresh_fail = quantile(exprs(chosen_units_eSet)[exprs(chosen_units_eSet) > 0], 0.2) #10
num_fail = 10

# Initial Gene Filtering
init.gf.vec = rowSums(exprs(chosen_units_eSet) > thresh_fail) > num_fail
chosen_units_eSet = chosen_units_eSet[init.gf.vec, ]
printf("Kept only %d genes expressed in more than %.2f units in more than %d cells , excluded %d genes\n", sum(init.gf.vec), thresh_fail, num_fail, sum(!init.gf.vec))






#The counts data frame contains feature-level read counts from tophat alignments of 96 single-cell libraries to the mm10 reference genome [@trapnell2009]. qc contains library alignment metrics obtained from Picard.
#de and hk are positive and negative control gene sets derived from population studies. batch and bio are factors labeling batch and time point respectively. Consider the joint distribution of these factors:
counts = exprs(chosen_units_eSet)
qc = pData(protocolData(chosen_units_eSet))

#prevent casing problems
rownames(counts) = sapply(rownames(counts), toupper)
housekeeping_genes = sapply(as.vector(as.matrix(read.table(housekeeping_list, stringsAsFactors = FALSE))), toupper)
hk = intersect(housekeeping_genes, rownames(counts))



#####################################
#start to do real work
#####################################
unnecessaryQC_metrics = apply(is.na(qc), 2, any) | apply(qc, 2, sd) < 1e-5
printf("remove %d QC metrics (%s) that either contained NAs or were constant across cells\n", sum(unnecessaryQC_metrics), paste(colnames(qc)[unnecessaryQC_metrics], collapse = "; "))
qc = qc[, !unnecessaryQC_metrics]

colnames(qc)

## Joint distribution of batches and biological conditions (time after induction)
if(!is.null(bio)) {
  table(batch,bio)
}

# Color scheme
cc <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"),brewer.pal(9,"Set3"))

# Barplot of read proportion mapping to the genome
barplot(qc$RALIGN[order(batch)], col=cc[batch][order(batch)], border=cc[batch][order(batch)], main="Percentage of mapped reads, colored by batch")
legend("bottomleft", legend=levels(batch), fill=cc,cex=0.4)

# Barplot of total read number
barplot(qc$NREADS[order(batch)], col=cc[batch][order(batch)], border=cc[batch][order(batch)], main="Total number of reads, colored by batch")
legend("topright", legend=levels(batch), fill=cc, cex=0.4)

qpc = prcomp(qc, center = TRUE, scale = TRUE)
barplot(cumsum(qpc$sdev^2)/sum(qpc$sdev^2), border="gray", xlab="PC", ylab="proportion of variance", main="Quality PCA")

barplot(qpc$x[,1][order(batch)], col=cc[batch][order(batch)], border=cc[batch][order(batch)], main="Quality PC1, colored by batch")
legend("bottomleft", legend=levels(batch), fill=cc, cex=0.8)


# Mean log10(x+1) expression
mu_obs = rowMeans(log10(counts[hk,]+1)) #0 if hk == NULL

# Drop-outs
drop_outs = counts[hk,] == 0

# Logistic Regression Model of Failure
logistic_coef = matrix(0,ncol(drop_outs),2)
for (si in seq_len(ncol(drop_outs))){
  fit = suppressWarnings(glm(cbind(drop_outs[,si],1 - drop_outs[,si]) ~ mu_obs,family=binomial(logit)))
  if(fit$converged){
    logistic_coef[si,] = fit$coefficients
  } else {
    logistic_coef[si,1] = logit(drop_rate[si])
  }
}


par(mfrow=c(1,2))

# Plot Failure Curves and Calculate AUC
plot(NULL, main = "False Negative Rate Curves",ylim = c(0,1),xlim = c(0,6), ylab = "Failure Probability", xlab = "Mean log10 Expression")
x = (0:60)/10
AUC = NULL
for(si in 1:ncol(counts)){
  y = 1/(exp(-logistic_coef[si,1] - logistic_coef[si,2] * x) + 1)
  AUC[si] = sum(y)/10
  lines(x, 1/(exp(-logistic_coef[si,1] - logistic_coef[si,2] * x) + 1), type = 'l', lwd = 2, col = cc[batch][si])
}

# Barplot of FNR AUC
barplot(AUC[order(batch)], col=cc[batch][order(batch)], border=cc[batch][order(batch)], main="FNR AUC, colored by batch")
legend("topleft", legend=levels(batch), fill=cc, cex=0.4)

mfilt_report <- metric_sample_filter(expr = counts,
                                     nreads = qc$NREADS,ralign = qc$RALIGN,
                                     suff_nreads = 10^5,
                                     suff_ralign = 90,
                                     pos_controls = rownames(counts) %in% hk,
                                     zcut = 3,mixture = FALSE, plot = TRUE)


hist(qc$RALIGN, breaks = 0:100)
# Hard threshold
abline(v = 15, col = "yellow", lwd = 2) 
# 3 (zcut) standard deviations below the mean ralign value
abline(v = mean(qc$RALIGN) - 3*sd(qc$RALIGN), col = "green", lwd = 2) 
# 3 (zcut) MADs below the median ralign value
abline(v = median(qc$RALIGN) - 3*mad(qc$RALIGN), col = "red", lwd = 2)
# Sufficient threshold
abline(v = 90, col = "grey", lwd = 2)
# Final threshold is the minimum of 1) the sufficient threshold and 2) the max of all others
thresh = min(90,max(c(15,mean(qc$RALIGN) - 3*sd(qc$RALIGN),median(qc$RALIGN) - 3*mad(qc$RALIGN))))
abline(v = thresh, col = "blue", lwd = 2, lty = 2)

legend("topleft",legend = c("Hard","SD","MAD","Sufficient","Final"),lwd = 2, col = c("yellow","green","red","grey","blue"),lty = c(1,1,1,1,2), cex = .5)


# Which thresholds are missing? (breadth)
is_na_filt = unlist(lapply(is.na(mfilt_report),any)) 

# Identify samples failing any threshold
m_sampfilter = !apply(simplify2array(mfilt_report[!is_na_filt]),1,any)

# Filter Samples
fcounts = counts[,m_sampfilter]
fqc = qc[m_sampfilter,]
fbatch = batch[m_sampfilter]
fbio = bio[m_sampfilter] #if bio==NULL then we get from this line fbio=NULL

# # Simple gene filter --> redundant with what I have above
# filterCount <- function(counts, nRead=5, nCell=5){
#   filter <- apply(counts, 1, function(x) length(x[x>=nRead])>=nCell)
#   return(filter)
# }
# genefilter <- filterCount(fcounts)
# fcounts = fcounts[genefilter,]

# skip filtering genes again
fcounts = counts
fhk = hk[hk %in% rownames(fcounts)]
fde = de[de %in% rownames(fcounts)] #if de==NULL then we get from this line fde=NULL


#FNR
#Michael: I see you recommend bulk_model=true even though the default is false? Why?
#I'm not clear on what this parameter does - what other gene features are possible if gFeatM is not specified?
fnr_out = estimate_ziber(x = fcounts, bulk_model = TRUE,
                         pos_controls = rownames(fcounts) %in% fhk,
                         verbose = TRUE, maxiter = 1000)


SUM_FN = function (ei) 
{
  sums = colSums(ei)
  eo = t(t(ei)*mean(sums)/sums)
  return(eo)
}

MED_FN_POS = function (ei) 
{
  zei = ei
  zei[ei == 0] = NA
  q = apply(zei, 2, quantile, 0.5, na.rm = TRUE)
  zeo = t(t(zei)/q) * mean(q)
  eo = zeo
  eo[ei == 0] = 0
  return(eo)
}

EFF_FN = function (ei) 
{
  sums = colSums(ei > 0)
  eo = t(t(ei)*sums/mean(sums))
  return(eo)
}

imputation = list(none=impute_null,expect=impute_expectation)
impute_args = list(p_nodrop = fnr_out$p_nodrop, mu = exp(fnr_out$Alpha[1,]))

scaling = list(none=identity,
             sum = SUM_FN,
             eff = EFF_FN,
             tmm = TMM_FN,
             uq = UQ_FN,
             uqp = UQ_FN_POS,
             fq = FQT_FN,
             fqp = FQ_FN_POS,
             deseq=DESEQ_FN,
             deseqp=DESEQ_FN_POS)



# Generate Params (RUN = FALSE)
#Michael: in your example qc=ppq which is ppq = scale(qc,center = TRUE,scale = TRUE)
#but I see in scone's code that it expects QC matrix and does the prcomp itself
params <- scone(expr = as.matrix(fcounts), scaling = scaling,
                imputation = imputation, impute_args = impute_args,
                ruv_negcon = fhk, k_ruv = 3,
                qc = as.matrix(fqc), k_qc = 3,
                bio = fbio, adjust_bio = adjust_bio,
                batch = fbatch, adjust_batch = adjust_batch,
                run = FALSE)
head(params)

apply(params,2,unique)

#exclide batch adjustment and no adjustment by bio and similar not-making-sense scenarios
is_screened = ((params$imputation_method == "none") & (params$scaling_method %in% c("deseq","uq","tmm"))) | 
  ((params$imputation_method == "expect") & (params$scaling_method %in% c("deseqp","uqp","fqp","eff"))) |
  ((params$adjust_biology == "bio") & (params$adjust_batch != "batch"))

params = params[!is_screened,]

##NO MORE PLOTS!
dev.off()

# Generate Scores and Ranking
tic = proc.time()
res = scone(expr = as.matrix(fcounts), scaling = scaling,
            imputation = imputation, impute_args = impute_args,
            ruv_negcon = fhk, k_ruv = 3,
            qc = as.matrix(fqc), k_qc = 3,
            bio = fbio, adjust_bio = adjust_bio,
            batch = fbatch, adjust_batch = adjust_batch,
            run = TRUE, params = params, verbose = TRUE,
            eval_poscon = fde, eval_kclust = 2:3, stratified_pam = stratified_pam,
            rezero = rezero, return_norm = "in_memory")
toc = proc.time()

printf("time elapsed:\n")
print(toc-tic)

save.image(file.path(out_dir, results_image_filename))


printf("Scone finished successfully and image saved")

}


#code to run this method and test it
if(FALSE) {
  rm(list=ls())
  #library(bioc2016singlecell)
  set.seed(112233) ## for reproducibility
  printf <- function(...) cat(sprintf(...))
  
  library(Biobase) #for AnnotatedDataFrame if scone is not loaded
  library(scone)
  library(RColorBrewer)
  
  BiocParallel::register(BiocParallel::SerialParam())
  
  setwd("/data/yosef2/Published_Data/TH17/code")
  source("scone_something.R")
  
  load("/data/yosef2/Published_Data/TH17/code/combined_th17_metadata.Rdata")
  housekeeping_list = "/data/yosef/CD8_effector_diff/src/SummaryPipeline/house_keeping_mouse_TitleCase.txt"
  out_dir = "/data/yosef2/Published_Data/TH17processed_aw20160712/scone/only_T16"
  
  #only T16 from the non-GFP mouse
  keep = combinedMetadata$michael_condition == "TGFB1_IL6-48h" &
    !combinedMetadata$is_population
  
  combinedMetadata = combinedMetadata[keep,]
  
  #tophaht + featureCounts
  collect_dir="/data/yosef2/Published_Data/TH17/processed_aw20160712/collect"
  load(file.path(collect_dir, "collectedRNASeqStudy.RData"))
  selected_eSet = collectedRNASeqStudy$cuff_eSet
  #featureNames(selected_eSet) = toupper(featureNames(selected_eSet))
  
  selected_eSet = selected_eSet[, combinedMetadata$SRR]
  batch = combinedMetadata$michael_batch
  scone_something(selected_eSet, housekeeping_list = housekeeping_list,
                  selected_expression_in_eSet = "counts_table",
                  batch = batch, bio = NULL, de = NULL,
                  adjust_bio = "no", adjust_batch = "yes", rezero = TRUE,
                  out_dir = getwd(), stratified_pam = TRUE,
                  results_image_filename = "sconeResults_tophatFeatureCounts.Rda",
                  plots_filename = "scone_plot_tophatFeatureCounts.pdf")
  
  stop("finished")
  
}


#code to inspect scone results
if(FALSE) {
  rm(list=ls())
  set.seed(4455)
  out_dir = "/home/eecs/allonwag/data2/TFH/processed_20160628/scone"
  load(file.path(out_dir, "sconeResults.Rda"))
  
  
  names(res)
  
  if(is.null(bio)) {
    #no positive control genes (de) so no evaluation of positive controls
    BIO_DEPENDENT_EVALUATION_METRICS = c("BIO_SIL", "EXP_WV_COR")
    res$scores = res$scores[, !(colnames(res$scores) %in% BIO_DEPENDENT_EVALUATION_METRICS)]
    res$metrics = res$metrics[, !(colnames(res$metrics) %in% BIO_DEPENDENT_EVALUATION_METRICS)]
  }
  
  head(res$scores)
  
  pc_obj = prcomp(res$scores[,-ncol(res$scores)],center = TRUE,scale = FALSE)
  bp_obj = biplot_colored(pc_obj,y = -res$scores[,ncol(res$scores)],expand = .6)
  
  #library(scone)
  #if I don't have this library installed
  source("biplot.R")
  bp_obj = biplot_colored(pc_obj,y = -res$scores[,ncol(res$scores)],expand = .6)
  
  points(bp_obj[grepl(",bio,batch",rownames(bp_obj)),], pch = 1, col = "red", cex = 1)
  points(bp_obj[grepl(",bio,batch",rownames(bp_obj)),], pch = 1, col = "red", cex = 1.5)
  
  
  bp_obj = biplot_colored(pc_obj,y = -res$scores[,ncol(res$scores)],expand = .6)
  
  points(bp_obj[grepl(",no_bio,batch",rownames(bp_obj)),], pch = 1, col = "red", cex = 1)
  points(bp_obj[grepl(",no_bio,batch",rownames(bp_obj)),], pch = 1, col = "red", cex = 1.5)
  
  #my addition to scone's tutorial
  bp_obj = biplot_colored(pc_obj,y = -res$scores[,ncol(res$scores)],expand = .6)
  
  points(bp_obj[grepl(",fq,",rownames(bp_obj)),], pch = 1, col = "green", cex = 1.5)
  points(bp_obj[grepl(",tmm,",rownames(bp_obj)),], pch = 1, col = "cyan", cex = 1.5)
  points(bp_obj[grepl(",deseq,",rownames(bp_obj)),], pch = 1, col = "magenta", cex = 1.5)
  points(bp_obj[grepl(",uqp,",rownames(bp_obj)),], pch = 1, col = "black", cex = 1.5)
  #end my addition
  
  bp_obj = biplot_colored(pc_obj,y = -res$scores[,ncol(res$scores)],expand = .6)
  
  points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1)
  points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1.5)
  
  points(t(bp_obj[rownames(bp_obj) == rownames(params)[1],]), pch = 1, col = "blue", cex = 1)
  points(t(bp_obj[rownames(bp_obj) == rownames(params)[1],]), pch = 1, col = "blue", cex = 1.5)
  
  arrows(bp_obj[rownames(bp_obj) == rownames(params)[1],][1],
         bp_obj[rownames(bp_obj) == rownames(params)[1],][2],
         bp_obj[1,][1],
         bp_obj[1,][2],
         lty = 2, lwd = 2)
  
  printf("(Trumpets): the selected method is %s", rownames(res$scores)[1])
  normalized_expression = res$normalized_data[[1]]
  selected_normalization_method = names(res$normalized_data)[1]
  save(normalized_expression, selected_normalization_method, file=file.path(out_dir, "sconeNormalizedExpression.Rda"))
  
  sessionInfo()
  
}
