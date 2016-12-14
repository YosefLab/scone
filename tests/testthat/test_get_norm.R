context("Tests for the get_normalized wrappers")
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

  # return_norm = hd5f
  res2 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=FALSE, run=TRUE, return_norm="hdf5", hdf5file = "tmp.h5")


  # return_norm = in_memory
  res3 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=FALSE, run=TRUE, return_norm="in_memory")

  #get_normalized can retrieve all norms
  all1 <- lapply(seq_along(assays(res3)), function(i) get_normalized(res3, i))
  all2 <- lapply(seq_along(assays(res3)), function(i) get_normalized(res3, names(assays(res3))[i]))

  expect_equal(all1, all2)

  #get_normalized should give the same results in three modes
  gn1 <- get_normalized(res1, 30)
  gn2 <- get_normalized(res2, 30)
  gn3 <- get_normalized(res3, 30)

  attributes(gn1) <- attributes(gn2) <- attributes(gn3) <- NULL
  expect_equal(gn1, gn2)
  expect_equal(gn1, gn3)

  gn1 <- get_normalized(res1, "none,uq,qc_k=2,bio,batch")
  gn2 <- get_normalized(res2, "none,uq,qc_k=2,bio,batch")
  gn3 <- get_normalized(res3, "none,uq,qc_k=2,bio,batch")

  attributes(gn1) <- attributes(gn2) <- attributes(gn3) <- NULL
  expect_equal(gn1, gn2)
  expect_equal(gn1, gn3)

  if (file.exists("tmp.h5")) file.remove("tmp.h5")
  }
)

test_that("get_normalized works in all three modes with nested model", {
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- as.factor(c(1,2,1,2,1,3,4,3,4,3))

  obj <- SconeExperiment(e, qc=qc_mat,
                         negcon_ruv=c(rep(TRUE, 100), rep(FALSE, NROW(e)-100)),
                         bio = as.factor(bio), batch=as.factor(batch))

  # return_norm = no
  res1 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=TRUE, eval_kclust = 2, run=TRUE)

  # return_norm = hd5f
  res2 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=TRUE, eval_kclust = 2, run=TRUE, return_norm="hdf5", hdf5file = "tmp.h5")


  # return_norm = in_memory
  res3 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=TRUE, eval_kclust = 2, run=TRUE, return_norm="in_memory")

  #get_normalized can retrieve all norms
  all1 <- lapply(seq_along(assays(res3)), function(i) get_normalized(res3, i))
  all2 <- lapply(seq_along(assays(res3)), function(i) get_normalized(res3, names(assays(res3))[i]))

  expect_equal(all1, all2)

  #get_normalized should give the same results in three modes
  gn1 <- get_normalized(res1, 30)
  gn2 <- get_normalized(res2, 30)
  gn3 <- get_normalized(res3, 30)

  attributes(gn1) <- attributes(gn2) <- attributes(gn3) <- NULL
  expect_equal(gn1, gn2)
  expect_equal(gn1, gn3)

  gn1 <- get_normalized(res1, "none,uq,qc_k=2,bio,batch")
  gn2 <- get_normalized(res2, "none,uq,qc_k=2,bio,batch")
  gn3 <- get_normalized(res3, "none,uq,qc_k=2,bio,batch")

  attributes(gn1) <- attributes(gn2) <- attributes(gn3) <- NULL
  expect_equal(gn1, gn2)
  expect_equal(gn1, gn3)

  expect_equal(get_scores(res1),get_scores(res2))
  expect_equal(get_scores(res1),get_scores(res3))

  if (file.exists("tmp.h5")) file.remove("tmp.h5")
}
)

test_that("get_normalized works with rezero", {
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
                evaluate=FALSE, run=TRUE, rezero=TRUE)

  # return_norm = in_memory
  res2 <- scone(obj, imputation=list(none=impute_null),
                scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                evaluate=FALSE, run=TRUE, return_norm="in_memory", rezero=TRUE)

  all1 <- lapply(seq_len(10), function(i) get_normalized(res1, i))
  all2 <- lapply(seq_len(10), function(i) get_normalized(res2, i))

  expect_equal(all1, all2)

}
)
