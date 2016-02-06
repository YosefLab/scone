source("../R/scone_main.R")
source("../R/zinb.R")
source("../R/SCONE_DEFAULTS.R")
source("../R/helper.R")
source("../R/scone_eval.R")
library(MASS)
library(RUVSeq)

e <-  matrix(rpois(10000, lambda = 5), ncol=10)
rownames(e) <- as.character(1:nrow(e))

# one combination
res <- scone(e, imputation=identity, scaling=identity, k_ruv=0, k_qc=0, evaluate=FALSE, run=TRUE)
stopifnot(all(res$normalized_data[[1]]==log1p(e)))

# add more imputations
res <- scone(e, imputation=list(a=identity, b=identity), scaling=identity, k_ruv=0, k_qc=0, evaluate=FALSE, run=TRUE)

# add more scaling
res <- scone(e, imputation=list(a=identity, b=identity), scaling=list(a=identity, b=identity, c=identity), k_ruv=0, k_qc=0, evaluate=FALSE, run=TRUE)

# add ruv (the first two should not work)
# scone(e, imputation=list(identity, identity), scaling=list(identity, identity, identity),
#       k_ruv=5, k_qc=0, evaluate=FALSE, run=FALSE)
# scone(e, imputation=list(identity, identity), scaling=list(identity, identity, identity),
#       k_ruv=5, ruv_negcon=1:100, k_qc=0, evaluate=FALSE, run=FALSE)
params <- scone(e, imputation=list(identity, identity), scaling=list(identity, identity, identity),
      k_ruv=5, ruv_negcon=as.character(1:100), k_qc=0, evaluate=FALSE, run=FALSE)

# add qc (the first two should not work)
qc_mat <- matrix(rnorm(20), nrow=10)
# scone(e, imputation=list(identity, identity), scaling=list(identity, identity, identity),
#       k_ruv=5, ruv_negcon=as.character(1:100), k_qc=5, evaluate=FALSE, run=FALSE)
# scone(e, imputation=list(identity, identity), scaling=list(identity, identity, identity),
#       k_ruv=5, ruv_negcon=as.character(1:100), k_qc=5, qc=qc_mat, evaluate=FALSE, run=FALSE)
res <- scone(e, imputation=list(identity, identity), scaling=list(identity, identity, identity),
      k_ruv=5, ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat, evaluate=FALSE, run=TRUE)

# add bio (the first two should not work)
bio <- rep(1:2, each=5)
# scone(e, imputation=list(identity, identity), scaling=list(identity, identity, identity),
#       k_ruv=5, ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
#       adjust_bio="yes", evaluate=FALSE, run=FALSE)
# scone(e, imputation=list(identity, identity), scaling=list(identity, identity, identity),
#       k_ruv=5, ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
#       adjust_bio="yes", bio=bio, evaluate=FALSE, run=FALSE)
bio <- as.factor(bio)
res <- scone(e, imputation=list(identity, identity), scaling=list(identity, identity, identity),
      k_ruv=5, ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
      adjust_bio="yes", bio=bio, evaluate=FALSE, run=TRUE)
res <- scone(e, imputation=list(identity, identity), scaling=list(identity, identity, identity),
      k_ruv=5, ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
      adjust_bio="force", bio=bio, evaluate=FALSE, run=TRUE)

# add batch (the first three should not work)
# batch <- rep(1:2, each=5)
# scone(e, imputation=list(identity, identity), scaling=list(identity, identity, identity),
#       k_ruv=5, ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
#       adjust_bio="force", bio=bio, adjust_batch="yes", evaluate=FALSE, run=FALSE)
# scone(e, imputation=list(identity, identity), scaling=list(identity, identity, identity),
#       k_ruv=5, ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
#       adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch, evaluate=FALSE, run=FALSE)
# batch <- as.factor(batch)
# scone(e, imputation=list(identity, identity), scaling=list(identity, identity, identity),
#       k_ruv=5, ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
#       adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch, evaluate=FALSE, run=FALSE)
batch <- as.factor(rep(1:2, 5))
res <- scone(e, imputation=list(a=identity, b=identity), scaling=list(a=identity, b=identity, c=identity),
      k_ruv=5, ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
      adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch, evaluate=FALSE, run=TRUE)

## real examples

# factorial
res <- scone(e, imputation=list(none=identity, zinb=impute_zinb), scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
             k_ruv=3, k_qc=2, ruv_negcon=as.character(1:100), qc=qc_mat,
             adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch, evaluate=FALSE, run=TRUE)

# nested
batch <- as.factor(c(1, 2, 1, 2, 1, 3, 4, 3, 4, 3))
params <- scone(e, imputation=list(none=identity, zinb=impute_zinb), scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
             k_ruv=3, k_qc=2, ruv_negcon=as.character(1:100), qc=qc_mat,
             adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch, evaluate=FALSE, run=FALSE)
params <- params[-(1:5),]
res <- scone(e, imputation=list(none=identity, zinb=impute_zinb), scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, ruv_negcon=as.character(1:100), qc=qc_mat,
                adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch, evaluate=FALSE, run=TRUE, params=params)

# evaluation
system.time(res <- scone(e, imputation=list(none=identity, zinb=impute_zinb), 
                         scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                         k_ruv=5, k_qc=2, ruv_negcon=as.character(1:100), qc=qc_mat,
                         adjust_bio="yes", bio=bio, adjust_batch="yes", batch=batch, run=TRUE,
                         evaluate=FALSE, eval_negcon=as.character(101:200), eval_poscon=as.character(201:300),
                         eval_knn=2, eval_kclust = 2, verbose=TRUE))

system.time(res <- scone(e, imputation=list(none=identity, zinb=impute_zinb), scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
             k_ruv=5, k_qc=2, ruv_negcon=as.character(1:100), qc=qc_mat,
             adjust_bio="yes", bio=bio, adjust_batch="yes", batch=batch, run=TRUE,
             evaluate=TRUE, eval_negcon=as.character(101:200), eval_poscon=as.character(201:300),
             eval_knn=2, eval_kclust = 2, verbose=TRUE))
