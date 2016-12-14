context("Tests for the get_design wrappers")
set.seed(421)
BiocParallel::register(BiocParallel::bpparam("SerialParam"))

test_that("get_normalized works in all three modes", {
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

  #get_normalized can retrieve all norms (slow)
  all1 <- lapply(seq_along(assays(res2)), function(i) get_design(res2, i))
  all2 <- lapply(seq_along(assays(res2)), function(i) get_design(res2, names(assays(res2))[i]))

  expect_equal(all1, all2)

  #get_normalized should give the same results in three modes
  gn1 <- get_design(res1, 30)
  gn2 <- get_design(res2, 30)

  expect_equal(gn1, gn2)

}
)
