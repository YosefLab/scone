context("Test the different back-ends of bplapply")
set.seed(13)

test_that("all back-ends work", {
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- as.factor(rep(1:2, 5))

  negcon_ruv <- c(rep(TRUE, 100), rep(FALSE, NROW(e)-100))

  obj <- sconeExperiment(e, bio=bio, batch=batch, qc=qc_mat, negcon_ruv=negcon_ruv)

  # serial
  res1 <- scone(obj, imputation=list(none=impute_null),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
               evaluate=TRUE, run=TRUE, return_norm = "in_memory",
               eval_kclust=2, bpparam=SerialParam())

  # multicore
  res2 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=TRUE, run=TRUE, return_norm = "in_memory",
                eval_kclust=2, bpparam=MulticoreParam(2))

  # snow
  res3 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=TRUE, run=TRUE, return_norm = "in_memory",
                eval_kclust=2, bpparam=SnowParam(workers=2, type="SOCK"))

  res4 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=TRUE, run=TRUE, return_norm = "in_memory",
                eval_kclust=2, bpparam=SnowParam(workers=2, type="FORK"))

  # batch jobs
  res5 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=TRUE, run=TRUE, return_norm = "in_memory",
                eval_kclust=2, bpparam=BatchJobsParam(2))

  if(require(doParallel)) {
    registerDoParallel(2)
    res6 <- scone(obj, imputation=list(none=impute_null),
                  scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                  k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                  evaluate=TRUE, run=TRUE, return_norm = "in_memory",
                  eval_kclust=2, bpparam=DoparParam())
    expect_equal(res1, res6)
  }

  expect_equal(res1, res2)
  expect_equal(res1, res3)
  expect_equal(res1, res4)
  expect_equal(res1, res5)
})


