context("Test the different back-ends of bplapply")
set.seed(13)

test_that("all back-ends work", {
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- as.factor(rep(1:2, 5))

  negcon_ruv <- c(rep(TRUE, 100), rep(FALSE, NROW(e)-100))

  obj <- SconeExperiment(e, bio=bio, batch=batch, qc=qc_mat, negcon_ruv=negcon_ruv)

  # serial
  res1 <- scone(obj, imputation=list(none=impute_null),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
               evaluate=TRUE, run=TRUE, return_norm = "in_memory",
               eval_kclust=2, bpparam=BiocParallel::SerialParam())

  # multicore
  if(.Platform$OS.type == "unix") {
    res2 <- scone(obj, imputation=list(none=impute_null),
                  scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                  k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                  evaluate=TRUE, run=TRUE, return_norm = "in_memory",
                  eval_kclust=2, bpparam=BiocParallel::MulticoreParam(2))

    res3 <- scone(obj, imputation=list(none=impute_null),
                  scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                  k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                  evaluate=TRUE, run=TRUE, return_norm = "in_memory",
                  eval_kclust=2, bpparam=BiocParallel::SnowParam(workers=2, type="FORK"))

    expect_equal(res1, res2)
    expect_equal(res1, res3)
  }

  # snow
  res4 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=TRUE, run=TRUE, return_norm = "in_memory",
                eval_kclust=2, bpparam=BiocParallel::SnowParam(workers=2, type="SOCK"))

  expect_equal(res1, res4)

  # batch jobs
  if(require(BatchJobs)) {
    res5 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=TRUE, run=TRUE, return_norm = "in_memory",
                eval_kclust=2, bpparam=BiocParallel::BatchJobsParam(2))
    expect_equal(res1, res5)
  }

  if(require(doParallel)) {
    registerDoParallel(2)
    res6 <- scone(obj, imputation=list(none=impute_null),
                  scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                  k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                  evaluate=TRUE, run=TRUE, return_norm = "in_memory",
                  eval_kclust=2, bpparam=BiocParallel::DoparParam())
    expect_equal(res1, res6)
  }

})


