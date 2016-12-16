context("Test normalization functions")
set.seed(1021)
BiocParallel::register(BiocParallel::bpparam("SerialParam"))

test_that("Upper-quartile normalization works the same as in the EDASeq package", {
  e <-  matrix(rpois(10000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))
  colnames(e) <- paste0("Sample", 1:ncol(e))

  negcon_ruv <- c(rep(TRUE, 100), rep(FALSE, NROW(e)-100))
  obj <- SconeExperiment(e, negcon_ruv=negcon_ruv)

  # UQ + RUV

  res <- scone(obj, imputation=impute_null, scaling=UQ_FN, k_ruv=5, k_qc=0,
               evaluate=FALSE, run=TRUE, return_norm = "in_memory")

  uq <- EDASeq::betweenLaneNormalization(e, which="upper", round=FALSE)
  rs <- lapply(1:5, function(i) RUVSeq::RUVg(log1p(uq), as.character(1:100), k=i, round=FALSE, isLog=TRUE)$norm)

  expect_equal(assay(res), uq)

  logres <- lapply(assays(res)[2:6], log1p)
  names(logres) <- NULL
  expect_equal(logres, rs)

  # UQ + QC
  qc_mat <- matrix(rnorm(20), nrow=10)
  obj <- SconeExperiment(e, negcon_ruv=negcon_ruv, qc=qc_mat)

  res <- scone(obj, imputation=impute_null, scaling=UQ_FN, k_ruv=0, k_qc=2,
               evaluate=FALSE, run=TRUE, return_norm = "in_memory")

  pca_qc <- prcomp(qc_mat, center=TRUE, scale=TRUE)
  qs <- lapply(1:2, function(i) {
    Y <- t(log1p(uq))
    W <- pca_qc$x[,1:i]
    alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
    correctedY <- Y - W %*% alpha
    return(t(correctedY))
  })

  logres <- lapply(assays(res)[2:3], log1p)
  names(logres) <- NULL
  expect_equal(logres, qs)
}
)
