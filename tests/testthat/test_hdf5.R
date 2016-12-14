context("Test the different return_norm options")
set.seed(13124323)
BiocParallel::register(BiocParallel::bpparam("SerialParam"))


test_that("hd5 checks", {
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- as.factor(rep(1:2, 5))
  negcon_ruv <- c(rep(TRUE, 100), rep(FALSE, NROW(e)-100))

  obj <- SconeExperiment(e, bio=bio, batch=batch, qc=qc_mat, negcon_ruv=negcon_ruv)

  # factorial
  expect_error(scone(obj, imputation=list(none=impute_null),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
               evaluate=FALSE, run=TRUE, return_norm = "hdf5"),
               "must be specified")

  expect_error(scone(obj, imputation=list(none=impute_null),
                     scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                     k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                     evaluate=FALSE, run=TRUE, return_norm = "hdf5", hdf5file = "/tmp/tmp", bpparam=BiocParallel::SnowParam(workers=2, type="SOCK")),
               "does not support multicores")
})

test_that("return_norm in memory", {
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- as.factor(rep(1:2, 5))

  negcon_ruv <- c(rep(TRUE, 100), rep(FALSE, NROW(e)-100))

  obj <- SconeExperiment(e, bio=bio, batch=batch, qc=qc_mat, negcon_ruv=negcon_ruv)

  # factorial
  res <- scone(obj, imputation=list(none=impute_null),
                     scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                     k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                     evaluate=FALSE, run=TRUE, return_norm = "in_memory")
  expect_true(length(assays(res))>1)
})

test_that("do not return_norm", {
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))
  colnames(e) <- paste0("Sample", 1:ncol(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- as.factor(rep(1:2, 5))

  negcon_ruv <- c(rep(TRUE, 100), rep(FALSE, NROW(e)-100))

  obj <- SconeExperiment(e, bio=bio, batch=batch, qc=qc_mat, negcon_ruv=negcon_ruv)

  # factorial
  res <- scone(obj, imputation=list(none=impute_null),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
               evaluate=FALSE, run=TRUE, return_norm = "no")
  expect_true(length(assays(res))==1)
  expect_equal(assays(res), assays(obj))
})
