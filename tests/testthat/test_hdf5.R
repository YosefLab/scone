context("Test the different return_norm options")
set.seed(13124323)
BiocParallel::register(bpparam("SerialParam"))


test_that("hd5 checks", {
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- as.factor(rep(1:2, 5))

  # factorial
  expect_error(scone(e, imputation=list(none=impute_null, zinb=impute_zinb),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=3, k_qc=2, ruv_negcon=as.character(1:100), qc=qc_mat,
               adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch,
               evaluate=FALSE, run=TRUE, return_norm = "hdf5"),
               "must be specified")

  BiocParallel::register(bpparam("MulticoreParam"))

  expect_error(scone(e, imputation=list(none=impute_null, zinb=impute_zinb),
                     scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                     k_ruv=3, k_qc=2, ruv_negcon=as.character(1:100), qc=qc_mat,
                     adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch,
                     evaluate=FALSE, run=TRUE, return_norm = "hdf5", hdf5file = "/tmp/tmp"),
               "does not support multicores")

  BiocParallel::register(bpparam("SerialParam"))
})

test_that("return_norm in memory", {
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- as.factor(rep(1:2, 5))

  # factorial
  res <- scone(e, imputation=list(none=impute_null, zinb=impute_zinb),
                     scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                     k_ruv=3, k_qc=2, ruv_negcon=as.character(1:100), qc=qc_mat,
                     adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch,
                     evaluate=FALSE, run=TRUE, return_norm = "in_memory")
  expect_is(res$normalized_data, "list")
  expect_true(all(sapply(res$normalized_data, is.matrix)))
})

test_that("do not return_norm", {
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- as.factor(rep(1:2, 5))

  # factorial
  res <- scone(e, imputation=list(none=impute_null, zinb=impute_zinb),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=3, k_qc=2, ruv_negcon=as.character(1:100), qc=qc_mat,
               adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch,
               evaluate=FALSE, run=TRUE, return_norm = "no")
  expect_null(res$normalized_data)
})
