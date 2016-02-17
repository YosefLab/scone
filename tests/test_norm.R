library(scone)
library(EDASeq)
library(RUVSeq)

e <-  matrix(rpois(10000, lambda = 5), ncol=10)
rownames(e) <- as.character(1:nrow(e))

# UQ + RUV
res <- scone(e, imputation=identity, scaling=UQ_FN, k_ruv=5, k_qc=0,
             evaluate=FALSE, run=TRUE, ruv_negcon=as.character(1:100))

uq <- betweenLaneNormalization(e, which="upper", round=FALSE)
rs <- lapply(1:5, function(i) RUVg(log1p(uq), as.character(1:100), k=i, round=FALSE, isLog=TRUE)$norm)
stopifnot(all(res$normalized_data[[2]]-rs[[1]]<1e-5))
stopifnot(all(res$normalized_data[[3]]-rs[[2]]<1e-5))
stopifnot(all(res$normalized_data[[4]]-rs[[3]]<1e-5))
stopifnot(all(res$normalized_data[[5]]-rs[[4]]<1e-5))
stopifnot(all(res$normalized_data[[6]]-rs[[5]]<1e-5))

# UQ + QC
qc_mat <- matrix(rnorm(20), nrow=10)

res <- scone(e, imputation=identity, scaling=UQ_FN, k_ruv=0, k_qc=2,
             evaluate=FALSE, run=TRUE, qc=qc_mat)

pca_qc <- prcomp(qc_mat, center=TRUE, scale=TRUE)
qs <- lapply(1:2, function(i) {
  Y <- t(log1p(uq))
  W <- pca_qc$x[,1:i]
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  correctedY <- Y - W %*% alpha
  return(t(correctedY))
  })

stopifnot(all(res$normalized_data[[2]]-qs[[1]]<1e-5))
stopifnot(all(res$normalized_data[[3]]-qs[[2]]<1e-5))
