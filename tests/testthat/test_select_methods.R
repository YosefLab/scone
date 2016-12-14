context("Tests for select_methods")
set.seed(421)
BiocParallel::register(BiocParallel::bpparam("SerialParam"))

test_that("select_methods works in all three modes", {
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- as.factor(rep(1:2, 5))

  obj <- SconeExperiment(e, qc=qc_mat,
                         negcon_ruv=c(rep(TRUE, 100), rep(FALSE, NROW(e)-100)),
                         bio = as.factor(bio), batch=as.factor(batch))

  # return_norm = no
  res1 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=FALSE, run=TRUE)

  # return_norm = in_memory
  res2 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=FALSE, run=TRUE, return_norm="in_memory")

  sub1 <- select_methods(res1, rownames(res1@scone_params)[1:5])
  sub2 <- select_methods(res1, 1:5)

  sub3 <- select_methods(res2, rownames(res2@scone_params)[1:5])
  sub4 <- select_methods(res2, 1:5)

  expect_equal(sub1, sub2)
  expect_equal(sub3, sub4)
  expect_equal(sub1@scone_params, sub3@scone_params)
})

test_that("get_normalized subsets score matrix", {
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- as.factor(rep(1:2, 5))

  obj <- SconeExperiment(e, qc=qc_mat,
                         negcon_ruv=c(rep(TRUE, 100), rep(FALSE, NROW(e)-100)),
                         bio = as.factor(bio), batch=as.factor(batch))

  # return_norm = no
  res1 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=TRUE, eval_kclust=2)

  sub1 <- select_methods(res1, rownames(res1@scone_params)[1:5])

  expect_equal(NROW(sub1@scone_scores), 5)
  expect_equal(NROW(sub1@scone_metrics), 5)

})
