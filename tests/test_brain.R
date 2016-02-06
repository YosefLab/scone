setwd("~/daviderisso/scone/tests")
source("../R/scone_main.R")
source("../R/zinb.R")
source("../R/SCONE_DEFAULTS.R")
source("../R/helper.R")
source("../R/scone_eval.R")
library(MASS)
library(RUVSeq)
library(brainUtils)
library(BiocParallel)

register(MulticoreParam(workers = 10))

### real data (remove from package)
input_dir <- "/data/yosef/BRAIN/processed_July2015/collect/pipe_davide/scone_input/"
counts <- as.matrix(read.table(file.path(input_dir, "filtered_counts.txt")))

set.seed(1999)
# idx_samples <- sample(1:ncol(counts), 200)
idx_samples <- 1:ncol(counts)
counts <- counts[, idx_samples]
filtered <- filterCount(counts, nRead = 10, nCell = 10)
qc <- as.matrix(read.table(file.path(input_dir, "qc.txt"))[idx_samples,])
batch <- read.table(file.path(input_dir, "batch.txt"))
batch <- droplevels(batch[idx_samples,1])
type <- read.table(file.path(input_dir, "cell_type.txt"))
type <- droplevels(type[idx_samples,1])

hk <- read.table("../data/house_keeping_mouse_TitleCase.txt")
hk <- intersect(rownames(filtered), hk[,1])

negcon <- read.table("../data/cell_cycle_Tirosh.txt", stringsAsFactors = FALSE)
nc <- rownames(filtered)[toupper(rownames(filtered)) %in% negcon[,2]]

poscon <- read.table("../data/macklis_markers.txt", stringsAsFactors = FALSE)
pc <- intersect(poscon[poscon[,2]!="Glia?",1], rownames(filtered))

params <- scone(filtered, imputation=list(none=identity),
             scaling=list(none=identity, fq=FQ_FN, deseq=DESEQ_FN), k_ruv=1, k_qc=1, 
             ruv_negcon=hk, qc=qc, adjust_bio="yes", bio=type, adjust_batch="yes", batch=batch,
             run=FALSE, evaluate=FALSE)
print(dim(params))
print(system.time(res <- scone(filtered, imputation=list(none=identity),
                         scaling=list(none=identity, fq=FQ_FN, deseq=DESEQ_FN), k_ruv=1, k_qc=1, 
                         ruv_negcon=hk, qc=qc, adjust_bio="yes", bio=type, adjust_batch="yes", batch=batch,
                         run=TRUE, evaluate=FALSE, verbose=TRUE)))
print(head(res$ranks))
rm(res)
print(system.time(res <- scone(filtered, imputation=list(none=identity),
                         scaling=list(none=identity, fq=FQ_FN, deseq=DESEQ_FN), k_ruv=1, k_qc=1, 
                         ruv_negcon=hk, qc=qc, adjust_bio="yes", bio=type, adjust_batch="yes", batch=batch,
                         run=TRUE, evaluate=TRUE, eval_negcon=nc, eval_poscon=pc,
                         verbose=TRUE
                         )))
print(head(res$ranks))
rm(res)
params <- scone(filtered, imputation=list(none=identity, zinb=impute_zinb),
                         scaling=list(none=identity, fq=FQ_FN, deseq=DESEQ_FN, tmm=TMM_FN), k_ruv=5, k_qc=5, 
                         ruv_negcon=hk, qc=qc, adjust_bio="yes", bio=type, adjust_batch="yes", batch=batch,
                         run=FALSE, evaluate=FALSE, eval_negcon=nc, eval_poscon=pc)
print(dim(params))

print(system.time(res <- scone(filtered, imputation=list(none=identity, zinb=impute_zinb),
                         scaling=list(none=identity, fq=FQ_FN, deseq=DESEQ_FN, tmm=TMM_FN), k_ruv=5, k_qc=5, 
                         ruv_negcon=hk, qc=qc, adjust_bio="yes", bio=type, adjust_batch="yes", batch=batch,
                         run=TRUE, evaluate=TRUE, eval_negcon=nc, eval_poscon=pc, verbose=TRUE
)))
print(head(res$ranks))
