context("General tests for the scone main function")
set.seed(13124)
BiocParallel::register(bpparam("SerialParam"))

test_that("Test with no real method (only identity)", {
  e <-  matrix(rpois(10000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  # one combination
  res <- scone(e, imputation=identity, scaling=identity, k_ruv=0, k_qc=0,
               evaluate=FALSE, run=TRUE)
  expect_equal(res$normalized_data[[1]], log1p(e))

  # add more imputations
  res <- scone(e, imputation=list(a=identity, b=identity), scaling=identity,
               k_ruv=0, k_qc=0, evaluate=FALSE, run=TRUE)
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
               evaluate=TRUE, run=TRUE, eval_knn=5, eval_kclust=5)

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
                           eval_knn=2, eval_kclust = 2, verbose=TRUE))

  expect_equal(rownames(res$evaluation), rownames(res$ranks))
  expect_equal(rownames(res$evaluation), rownames(res$params))
  expect_equal(rownames(res$evaluation), names(res$normalized_data))

  res2 <- scone(e, imputation=list(none=identity, zinb=impute_zinb),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=5, k_qc=2, ruv_negcon=as.character(1:10),
               qc=qc_mat, adjust_bio="yes", bio=bio,
               adjust_batch="yes", batch=batch, run=TRUE,
               evaluate=FALSE, eval_negcon=as.character(11:20),
               eval_poscon=as.character(21:30),
               eval_knn=2, eval_kclust = 2, verbose=TRUE)

  norm_ordered <- res2$normalized_data[names(res$normalized_data)]
  expect_equal(norm_ordered, res$normalized_data)
  expect_null(res2$evaluation)
  expect_null(res2$ranks)
})

test_that("scone works with only one normalization",{
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))

  res <- scone(e, imputation=list(none=identity),
               scaling=list(none=identity),
               k_ruv=0, k_qc=0, run=TRUE,
               evaluate=TRUE, eval_knn=2, eval_kclust = 2)

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
               k_ruv=0, k_qc=0, adjust_bio="yes", bio=bio,
               adjust_batch="yes", batch=batch, run=FALSE,
               evaluate=TRUE, eval_negcon=as.character(11:20),
               eval_poscon=as.character(21:30),
               eval_knn=2, eval_kclust = 2, conditional_pam = TRUE)

  expect_error(res <- scone(e, imputation=list(none=identity, zinb=impute_zinb),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=0, k_qc=0, adjust_bio="yes", bio=bio,
               adjust_batch="yes", batch=batch, run=FALSE,
               evaluate=TRUE, eval_negcon=as.character(11:20),
               eval_poscon=as.character(21:30),
               eval_knn=2, eval_kclust = 6, conditional_pam = TRUE),
               "For conditional_pam, max 'eval_kclust' must be smaller than min bio class size")

  expect_error(res <- scone(e, imputation=list(none=identity, zinb=impute_zinb),
                            scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                            k_ruv=0, k_qc=0,
                            adjust_batch="yes", batch=batch, run=FALSE,
                            evaluate=TRUE, eval_negcon=as.character(11:20),
                            eval_poscon=as.character(21:30),
                            eval_knn=2, eval_kclust = 6, conditional_pam = TRUE),
               "If `bio` is null, `conditional_pam` cannot be TRUE")

})
