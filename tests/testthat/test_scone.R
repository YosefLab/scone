context("General tests for the scone main function")
set.seed(13124)
BiocParallel::register(bpparam("SerialParam"))

test_that("Test with no real method (only identity)", {
  e <-  matrix(rpois(10000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  # one combination
  res <- scone(e, imputation=identity, scaling=identity, k_ruv=0, k_qc=0,
               evaluate=FALSE, run=TRUE, return_norm = "in_memory")
  expect_equal(res$normalized_data[[1]], log1p(e))

  # add more imputations
  res <- scone(e, imputation=list(a=identity, b=identity), scaling=identity,
               k_ruv=0, k_qc=0, evaluate=FALSE, run=TRUE, return_norm = "in_memory")
  expect_equal(res$normalized_data[[1]], log1p(e))
  expect_equal(res$normalized_data[[2]], log1p(e))

  # add more scaling
  res <- scone(e, imputation=list(a=identity, b=identity),
               scaling=list(a=identity, b=identity, c=identity), k_ruv=0,
               k_qc=0, evaluate=FALSE, run=TRUE)

  # add ruv (the first two should not work)
  expect_error(scone(e, imputation=list(identity, identity),
                     scaling=list(identity, identity, identity),
                     k_ruv=5, k_qc=0, evaluate=FALSE, run=FALSE),
        "ruv_negcon must be specified")

  expect_error(scone(e, imputation=list(identity, identity),
                     scaling=list(identity, identity, identity), k_ruv=5,
                     ruv_negcon=1:100, k_qc=0, evaluate=FALSE, run=FALSE),
        "must be a character vector")

  params <- scone(e, imputation=list(identity, identity),
                  scaling=list(identity, identity, identity), k_ruv=5,
                  ruv_negcon=as.character(1:100), k_qc=0, evaluate=FALSE,
                  run=FALSE)

  # add qc (the first two should not work)
  qc_mat <- matrix(rnorm(20), nrow=10)

  expect_error(scone(e, imputation=list(identity, identity),
        scaling=list(identity, identity, identity), k_ruv=5,
        ruv_negcon=as.character(1:100), k_qc=5, evaluate=FALSE, run=FALSE),
        "qc must be specified")
  expect_error(scone(e, imputation=list(identity, identity),
        scaling=list(identity, identity, identity), k_ruv=5,
        ruv_negcon=as.character(1:100), k_qc=5, qc=qc_mat, evaluate=FALSE,
        run=FALSE),
        "k_qc must be less or equal than the number of columns")

  res <- scone(e, imputation=list(identity, identity),
               scaling=list(identity, identity, identity), k_ruv=5,
               ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
               evaluate=FALSE, run=TRUE)

  # add bio (the first two should not work)
  bio <- rep(1:2, each=5)

  expect_error(scone(e, imputation=list(identity, identity),
                     scaling=list(identity, identity, identity), k_ruv=5,
                     ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
                     adjust_bio="yes", evaluate=FALSE, run=FALSE),
               "if adjust_bio is 'yes' or 'force', 'bio' must be specified")

  expect_error(scone(e, imputation=list(identity, identity),
                     scaling=list(identity, identity, identity), k_ruv=5,
                     ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
                     adjust_bio="yes", bio=bio, evaluate=FALSE, run=FALSE),
               "'bio' must be a factor")

  bio <- as.factor(bio)
  res <- scone(e, imputation=list(identity, identity),
               scaling=list(identity, identity, identity), k_ruv=5,
               ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
               adjust_bio="yes", bio=bio, evaluate=FALSE, run=TRUE)
  res <- scone(e, imputation=list(identity, identity),
               scaling=list(identity, identity, identity), k_ruv=5,
               ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
               adjust_bio="force", bio=bio, evaluate=FALSE, run=TRUE)

  # add batch (the first three should not work)
  batch <- rep(1:2, each=5)

  expect_error(scone(e, imputation=list(identity, identity),
                     scaling=list(identity, identity, identity), k_ruv=5,
                     ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
                     adjust_bio="force", bio=bio, adjust_batch="yes",
                     evaluate=FALSE, run=FALSE),
               "if adjust_batch is 'yes' or 'force', 'batch' must be specified")

  expect_error(scone(e, imputation=list(identity, identity),
                     scaling=list(identity, identity, identity), k_ruv=5,
                     ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
                     adjust_bio="force", bio=bio, adjust_batch="yes",
                     batch=batch, evaluate=FALSE, run=FALSE),
               "'batch' must be a factor")

  batch <- as.factor(batch)
  expect_error(scone(e, imputation=list(identity, identity),
                     scaling=list(identity, identity, identity), k_ruv=5,
                     ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
                     adjust_bio="force", bio=bio, adjust_batch="yes",
                     batch=batch, evaluate=FALSE, run=FALSE),
               "Biological conditions and batches are confounded")

  batch <- as.factor(rep(1:2, 5))
  res <- scone(e, imputation=list(a=identity, b=identity),
               scaling=list(a=identity, b=identity, c=identity),
               k_ruv=5, ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
               adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch,
               evaluate=FALSE, run=TRUE)

  ## add evaluation
  res <- scone(e, imputation=list(a=identity, b=identity),
               scaling=list(a=identity, b=identity, c=identity),
               k_ruv=5, ruv_negcon=as.character(1:100), k_qc=2, qc=qc_mat,
               adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch,
               evaluate=TRUE, run=TRUE, eval_kclust=5)

})

test_that("Test imputation and scaling", {
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- as.factor(rep(1:2, 5))

  # factorial
  res <- scone(e, imputation=list(none=identity, zinb=impute_zinb),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=3, k_qc=2, ruv_negcon=as.character(1:100), qc=qc_mat,
               adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch,
               evaluate=FALSE, run=TRUE)

  # nested
  batch <- as.factor(c(1, 2, 1, 2, 1, 3, 4, 3, 4, 3))
  params <- scone(e, imputation=list(none=identity, zinb=impute_zinb),
                  scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                  k_ruv=3, k_qc=2, ruv_negcon=as.character(1:10), qc=qc_mat,
                  adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch,
                  evaluate=FALSE, run=FALSE)
  params <- params[-(1:5),]
  res <- scone(e, imputation=list(none=identity, zinb=impute_zinb),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=3, k_qc=2, ruv_negcon=as.character(1:10), qc=qc_mat,
               adjust_bio="force", bio=bio, adjust_batch="yes", batch=batch,
               evaluate=FALSE, run=TRUE, params=params)

  # evaluation
  system.time(res <- scone(e, imputation=list(none=identity, zinb=impute_zinb),
                           scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                           k_ruv=5, k_qc=2, ruv_negcon=as.character(1:10),
                           qc=qc_mat, adjust_bio="yes", bio=bio,
                           adjust_batch="yes", batch=batch, run=TRUE,
                           evaluate=TRUE, eval_negcon=as.character(11:20),
                           eval_poscon=as.character(21:30),
                           eval_kclust = 2, verbose=FALSE,
                           return_norm = "in_memory"))

  expect_equal(rownames(res$metrics), rownames(res$scores))
  expect_equal(rownames(res$metrics), rownames(res$params))
  expect_equal(rownames(res$metrics), names(res$normalized_data))

  res2 <- scone(e, imputation=list(none=identity, zinb=impute_zinb),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=5, k_qc=2, ruv_negcon=as.character(1:10),
               qc=qc_mat, adjust_bio="yes", bio=bio,
               adjust_batch="yes", batch=batch, run=TRUE,
               evaluate=FALSE, eval_negcon=as.character(11:20),
               eval_poscon=as.character(21:30),
               eval_kclust = 2, verbose=FALSE, return_norm = "in_memory")

  norm_ordered <- res2$normalized_data[names(res$normalized_data)]
  expect_equal(norm_ordered, res$normalized_data)
  expect_null(res2$metrics)
  expect_null(res2$scores)
})

test_that("scone works with only one normalization",{
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  res <- scone(e, imputation=list(none=identity),
               scaling=list(none=identity),
               k_ruv=0, k_qc=0, run=TRUE,
               evaluate=TRUE, eval_kclust = 2, return_norm = "in_memory")

  expect_equal(res$normalized_data[[1]], log1p(e))
})

test_that("conditional PAM",{
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- gl(5, 2)

  res <- scone(e, imputation=list(none=identity, zinb=impute_zinb),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=0, k_qc=0, adjust_bio="yes", bio=bio, run=FALSE,
               evaluate=TRUE, eval_negcon=as.character(11:20),
               eval_poscon=as.character(21:30),
               eval_kclust = 2, stratified_pam = TRUE)

  expect_error(res <- scone(e, imputation=list(none=identity, zinb=impute_zinb),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=0, k_qc=0, adjust_bio="yes", bio=bio,
               adjust_batch="yes", batch=batch, run=FALSE,
               evaluate=TRUE, eval_negcon=as.character(11:20),
               eval_poscon=as.character(21:30),
               eval_kclust = 6, stratified_pam = TRUE),
               "For stratified_pam, max 'eval_kclust' must be smaller than bio-cross-batch stratum size")

  expect_error(res <- scone(e, imputation=list(none=identity, zinb=impute_zinb),
                            scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                            k_ruv=0, k_qc=0, run = FALSE,
                            evaluate = TRUE, eval_negcon=as.character(11:20),
                            eval_poscon=as.character(21:30),
                            eval_kclust = 6, stratified_pam = TRUE),
               "For stratified_pam, bio and/or batch must be specified")

})


test_that("if bio=no bio is ignored", {
  e <-  matrix(rpois(10000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  res1 <- scone(e, imputation=identity, scaling=identity, k_ruv=0, k_qc=0,
               adjust_bio = "no", bio=gl(2, 5), eval_kclust = 3,
               return_norm = "in_memory")

  res2 <- scone(e, imputation=identity, scaling=identity, k_ruv=0, k_qc=0,
                 adjust_bio = "no", eval_kclust = 3,
                return_norm = "in_memory")

  expect_equal(res1$normalized_data, res2$normalized_data)
})

test_that("if batch=no batch is ignored", {
  e <-  matrix(rpois(10000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  res1 <- scone(e, imputation=identity, scaling=identity, k_ruv=0, k_qc=0,
                adjust_batch = "no", batch=gl(2, 5), eval_kclust = 3,
                return_norm = "in_memory")

  res2 <- scone(e, imputation=identity, scaling=identity, k_ruv=0, k_qc=0,
                adjust_batch = "no", eval_kclust = 3,
                return_norm = "in_memory")

  expect_equal(res1$normalized_data, res2$normalized_data)
})

test_that("batch and bio can be confounded if at least one of adjust_bio or adjust_batch is no", {
  e <-  matrix(rpois(10000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  expect_warning(scone(e, imputation=identity, scaling=identity, k_ruv=0, k_qc=0,
                adjust_batch = "no", batch=gl(2, 5), bio=gl(2, 5), eval_kclust = 3),
                "Biological conditions and batches are confounded.")

  expect_warning(scone(e, imputation=identity, scaling=identity, k_ruv=0, k_qc=0,
                       adjust_bio = "no", batch=gl(2, 5), bio=gl(2, 5), eval_kclust = 3),
                 "Biological conditions and batches are confounded.")
})

test_that("batch and bio can contain NA", {
  e <-  matrix(rpois(10000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  batch <- gl(2, 5)
  bio <- gl(5, 2)

  res1 <- scone(e, imputation=identity, scaling=identity, k_ruv=0, k_qc=0, evaluate = TRUE,
               adjust_batch = "no", batch=batch, bio=bio, eval_kclust = 3)

  batch[1] <- NA
  bio[2] <- NA

  res2 <- scone(e, imputation=identity, scaling=identity, k_ruv=0, k_qc=0, evaluate = TRUE,
        adjust_batch = "no", batch=batch, bio=bio, eval_kclust = 3)

  expect_true(!is.na(res2$metrics[,"BIO_SIL"]))
  expect_true(!is.na(res2$metrics[,"BATCH_SIL"]))
  expect_false(all(res1$metrics[,"BIO_SIL"] == res2$metrics[,"BIO_SIL"]))
  expect_false(all(res1$metrics[,"BATCH_SIL"] == res2$metrics[,"BATCH_SIL"]))
})
