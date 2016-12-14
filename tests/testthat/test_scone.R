context("General tests for the scone main function")
set.seed(13124)
BiocParallel::register(BiocParallel::bpparam("SerialParam"))

test_that("Test with no real method (only identity)", {
  e <-  matrix(rpois(10000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))
  colnames(e) <- paste0("Sample", 1:ncol(e))
  obj <- SconeExperiment(e)

  # one combination
  res <- scone(obj, imputation=impute_null, scaling=identity, k_ruv=0, k_qc=0,
               evaluate=FALSE, run=TRUE, return_norm = "in_memory")
  expect_equal(assay(res), e)

  # res2 should be the same as res
  res2 <- scone(res, imputation=impute_null, scaling=identity, k_ruv=0, k_qc=0,
               evaluate=FALSE, run=TRUE, return_norm = "in_memory")
  expect_equal(res2, res)

  # add more imputations
  res <- scone(obj, imputation=list(a=impute_null, b=impute_null), scaling=identity,
               k_ruv=0, k_qc=0, evaluate=FALSE, run=TRUE, return_norm = "in_memory")
  expect_equal(assay(res), e)
  expect_equal(assay(res, 2), e)

  # add more scaling
  res <- scone(obj, imputation=list(a=impute_null, b=impute_null),
               scaling=list(a=identity, b=identity, c=identity), k_ruv=0,
               k_qc=0, evaluate=FALSE, run=TRUE)

  # add ruv
  expect_error(scone(obj, imputation=list(impute_null,impute_null),
                     scaling=list(identity, identity, identity),
                     k_ruv=5, k_qc=0, evaluate=FALSE, run=FALSE),
               "negative controls must be specified")

  obj <- SconeExperiment(e, negcon_ruv=c(rep(TRUE, 100), rep(FALSE, NROW(e)-100)))
  obj <- scone(obj, imputation=list(impute_null,impute_null),
               scaling=list(identity, identity, identity), k_ruv=5,
               k_qc=0, evaluate=FALSE, run=FALSE)
  obj2 <- scone(obj, imputation=list(impute_null,impute_null),
                scaling=list(identity, identity, identity), k_ruv=5,
                k_qc=0, evaluate=FALSE, run=FALSE)
  expect_equal(obj, obj2)

  # add qc
  expect_error(scone(obj, imputation=list(impute_null,impute_null),
                     scaling=list(identity, identity, identity), k_ruv=5,
                     k_qc=5, evaluate=FALSE, run=FALSE),
               "QC metrics must be specified")

  qc_mat <- matrix(rnorm(20), nrow=10)
  obj <- SconeExperiment(e, qc=qc_mat, negcon_ruv=c(rep(TRUE, 100), rep(FALSE, NROW(e)-100)))

  res <- scone(obj, imputation=list(impute_null,impute_null),
               scaling=list(identity, identity, identity), k_ruv=5, k_qc=2,
               evaluate=FALSE, run=TRUE)

  # add bio
  bio <- rep(1:2, each=5)

  expect_error(scone(obj, imputation=list(impute_null,impute_null),
                     scaling=list(identity, identity, identity), k_ruv=5,
                     k_qc=2, adjust_bio="yes", evaluate=FALSE, run=FALSE),
               "if adjust_bio is 'yes' or 'force', 'bio' must be specified")

  obj <- SconeExperiment(e, qc=qc_mat,
                         negcon_ruv=c(rep(TRUE, 100), rep(FALSE, NROW(e)-100)),
                         bio = as.factor(bio))

  res <- scone(obj, imputation=list(impute_null,impute_null),
               scaling=list(identity, identity, identity), k_ruv=5,
               k_qc=2, adjust_bio="yes", evaluate=FALSE, run=TRUE)
  res <- scone(obj, imputation=list(impute_null,impute_null),
               scaling=list(identity, identity, identity), k_ruv=5,
               k_qc=2, adjust_bio="force", evaluate=FALSE, run=TRUE)

  # add batch
  batch <- rep(1:2, each=5)

  expect_error(scone(obj, imputation=list(impute_null,impute_null),
                     scaling=list(identity, identity, identity), k_ruv=5,
                     k_qc=2, adjust_bio="force", adjust_batch="yes",
                     evaluate=FALSE, run=FALSE),
               "if adjust_batch is 'yes' or 'force', 'batch' must be specified")

  obj <- SconeExperiment(e, qc=qc_mat,
                         negcon_ruv=c(rep(TRUE, 100), rep(FALSE, NROW(e)-100)),
                         bio = as.factor(bio), batch=as.factor(batch))

  expect_error(scone(obj, imputation=list(impute_null,impute_null),
                     scaling=list(identity, identity, identity), k_ruv=5,
                     k_qc=2, adjust_bio="force", adjust_batch="yes",
                     evaluate=FALSE, run=FALSE),
               "Biological conditions and batches are confounded")

  batch <- as.factor(rep(1:2, 5))
  obj <- SconeExperiment(e, qc=qc_mat,
                         negcon_ruv=c(rep(TRUE, 100), rep(FALSE, NROW(e)-100)),
                         bio = as.factor(bio), batch=as.factor(batch))
  res <- scone(obj, imputation=list(a=impute_null, b=impute_null),
               scaling=list(a=identity, b=identity, c=identity),
               k_ruv=5, k_qc=2, adjust_bio="force", adjust_batch="yes",
               evaluate=FALSE, run=TRUE)

  ## add evaluation
  res <- scone(obj, imputation=list(a=impute_null, b=impute_null),
               scaling=list(a=identity, b=identity, c=identity),
               k_ruv=5, k_qc=2, adjust_bio="force", adjust_batch="yes",
               evaluate=TRUE, run=TRUE, eval_kclust=5)

})

test_that("Test imputation and scaling", {
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))
  colnames(e) <- paste0("Sample", 1:ncol(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- as.factor(rep(1:2, 5))

  obj <- SconeExperiment(e, qc=qc_mat,
                         negcon_ruv=c(rep(TRUE, 100), rep(FALSE, NROW(e)-100)),
                         bio = as.factor(bio), batch=as.factor(batch))

  # factorial
  res <- scone(obj, imputation=list(none=impute_null),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
               evaluate=FALSE, run=TRUE)

  # nested
  batch <- as.factor(c(1, 2, 1, 2, 1, 3, 4, 3, 4, 3))
  obj <- SconeExperiment(e, qc=qc_mat,
                         negcon_ruv=c(rep(TRUE, 100), rep(FALSE, NROW(e)-100)),
                         bio = as.factor(bio), batch=as.factor(batch))

  obj <- scone(obj, imputation=list(none=impute_null),
                  scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                  k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
                  evaluate=FALSE, run=FALSE)
  obj@scone_params <- obj@scone_params[-(1:5),]

  res <- scone(obj, imputation=list(none=impute_null),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=3, k_qc=2, adjust_bio="force", adjust_batch="yes",
               evaluate=FALSE, run=TRUE)

  # evaluation
  ruv_negcon <- eval_negcon <- eval_poscon <- rep(FALSE, NROW(e))
  ruv_negcon[1:10] <- TRUE
  eval_negcon[11:20] <- TRUE
  eval_poscon[21:30] <- TRUE
  obj <- SconeExperiment(e, qc=qc_mat, negcon_ruv=ruv_negcon,
                         negcon_eval=eval_negcon, poscon=eval_poscon,
                         bio=as.factor(bio), batch=as.factor(batch))

  system.time(res <- scone(obj, imputation=list(none=impute_null),
                           scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                           k_ruv=5, k_qc=2, adjust_bio="yes",
                           adjust_batch="yes", run=TRUE, evaluate=TRUE,
                           eval_kclust = 2, verbose=FALSE,
                           return_norm = "in_memory"))

  expect_equal(rownames(res@scone_metrics), rownames(res@scone_scores))
  expect_equal(rownames(res@scone_metrics), rownames(res@scone_params))
  expect_equal(rownames(res@scone_metrics), names(assays(res)))

  res2 <- scone(obj, imputation=list(none=impute_null),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=5, k_qc=2, adjust_bio="yes", adjust_batch="yes", run=TRUE,
               evaluate=FALSE, eval_kclust = 2, verbose=FALSE,
               return_norm = "in_memory")

  norm_ordered <- assays(res2)[names(assays(res))]
  expect_equal(norm_ordered, assays(res))
  expect_true(is.na(res2@scone_scores))
  expect_true(is.na(res2@scone_metrics))
})

test_that("scone works with only one normalization",{
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))
  colnames(e) <- paste0("Sample", 1:ncol(e))
  obj <- SconeExperiment(e)

  res <- scone(obj, imputation=list(none=impute_null),
               scaling=list(none=identity),
               k_ruv=0, k_qc=0, run=TRUE,
               evaluate=TRUE, eval_kclust = 2, return_norm = "in_memory")

  expect_equal(assay(res), e)
})

test_that("conditional PAM",{
  e <-  matrix(rpois(1000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))
  colnames(e) <- paste0("Sample", 1:ncol(e))

  qc_mat <- matrix(rnorm(20), nrow=10)
  bio <- gl(2, 5)
  batch <- gl(5, 2)
  eval_negcon <- eval_poscon <- rep(FALSE, NROW(e))
  eval_negcon[11:20] <- TRUE
  eval_poscon[21:30] <- TRUE

  obj <- SconeExperiment(e, qc=qc_mat, bio=bio,
                         negcon_eval = eval_negcon, poscon=eval_poscon)

  res <- scone(obj, imputation=list(none=impute_null),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=0, k_qc=0, adjust_bio="yes", run=FALSE,
               evaluate=TRUE, eval_kclust = 2, stratified_pam = TRUE)

  obj <- SconeExperiment(e, qc=qc_mat, bio=bio, batch=batch,
                         negcon_eval = eval_negcon, poscon=eval_poscon)

  expect_error(res <- scone(obj, imputation=list(none=impute_null),
               scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
               k_ruv=0, k_qc=0, adjust_bio="yes", adjust_batch="yes", run=FALSE,
               evaluate=TRUE, eval_kclust = 6, stratified_pam = TRUE),
               "For stratified_pam, max 'eval_kclust' must be smaller than bio-cross-batch stratum size")

  obj <- SconeExperiment(e, qc=qc_mat, negcon_eval = eval_negcon, poscon=eval_poscon)

  expect_error(res <- scone(obj, imputation=list(none=impute_null),
                            scaling=list(none=identity, uq=UQ_FN, deseq=DESEQ_FN),
                            k_ruv=0, k_qc=0, run = FALSE,
                            evaluate = TRUE,
                            eval_kclust = 6, stratified_pam = TRUE),
               "For stratified_pam, bio and/or batch must be specified")

})


test_that("if bio=no bio is ignored", {
  e <-  matrix(rpois(10000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))
  colnames(e) <- paste0("Sample", 1:ncol(e))
  bio <- gl(2, 5)
  obj1 <- SconeExperiment(e)
  obj2 <- SconeExperiment(e, bio=bio)

  res1 <- scone(obj1, imputation=impute_null, scaling=identity, k_ruv=0, k_qc=0,
               adjust_bio = "no",  eval_kclust = 3, return_norm = "in_memory")

  res2 <- scone(obj2, imputation=impute_null, scaling=identity, k_ruv=0, k_qc=0,
                 adjust_bio = "no", eval_kclust = 3, return_norm = "in_memory")

  expect_equal(assays(res1), assays(res2))
})

test_that("if batch=no batch is ignored", {
  e <-  matrix(rpois(10000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))
  colnames(e) <- paste0("Sample", 1:ncol(e))
  batch <- gl(2, 5)
  obj1 <- SconeExperiment(e)
  obj2 <- SconeExperiment(e, batch=batch)

  res1 <- scone(obj1, imputation=impute_null, scaling=identity, k_ruv=0, k_qc=0,
                adjust_batch = "no", eval_kclust = 3, return_norm = "in_memory")

  res2 <- scone(obj2, imputation=impute_null, scaling=identity, k_ruv=0, k_qc=0,
                adjust_batch = "no", eval_kclust = 3, return_norm = "in_memory")

  expect_equal(assays(res1), assays(res2))
})

test_that("batch and bio can be confounded if at least one of adjust_bio or adjust_batch is no", {
  e <-  matrix(rpois(10000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))
  colnames(e) <- paste0("Sample", 1:ncol(e))
  obj <- SconeExperiment(e, batch=gl(2, 5), bio=gl(2, 5))

  expect_warning(scone(obj, imputation=impute_null, scaling=identity, k_ruv=0, k_qc=0,
                adjust_batch = "yes", eval_kclust = 3),
                "Biological conditions and batches are confounded.")

  expect_warning(scone(obj, imputation=impute_null, scaling=identity, k_ruv=0, k_qc=0,
                       adjust_bio = "yes", eval_kclust = 3),
                 "Biological conditions and batches are confounded.")
})

test_that("batch and bio can contain NA", {
  e <-  matrix(rpois(10000, lambda = 5), ncol=10)
  rownames(e) <- as.character(1:nrow(e))
  colnames(e) <- paste0("Sample", 1:ncol(e))
  batch <- gl(2, 5)
  bio <- gl(5, 2)
  obj <- SconeExperiment(e, batch=batch, bio=bio)
  res1 <- scone(obj, imputation=impute_null, scaling=identity, k_ruv=0, k_qc=0, evaluate = TRUE,
               adjust_batch = "no", eval_kclust = 3)

  batch[1] <- NA
  bio[2] <- NA
  obj <- SconeExperiment(e, batch=batch, bio=bio)

  res2 <- scone(obj, imputation=impute_null, scaling=identity, k_ruv=0, k_qc=0, evaluate = TRUE,
        adjust_batch = "no", eval_kclust = 3)

  expect_true(!is.na(res2@scone_metrics[,"BIO_SIL"]))
  expect_true(!is.na(res2@scone_metrics[,"BATCH_SIL"]))
  expect_false(all(res1@scone_metrics[,"BIO_SIL"] == res2@scone_metrics[,"BIO_SIL"]))
  expect_false(all(res1@scone_metrics[,"BATCH_SIL"] == res2@scone_metrics[,"BATCH_SIL"]))
})
