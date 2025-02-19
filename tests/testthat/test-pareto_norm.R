test_that("PSiNorm works with all input classes", {
  m <- matrix(0L, nrow=4, ncol=5)
  m[c(1,3) , 1] <- 1:2
  m[1:3, 3] <- c(2L, 9L, 3L)
  m[ , 4] <- 2:1
  m[2, 5] <- 1L

  out1 <- PsiNorm(m)

  svt <- as(m, "SVT_SparseMatrix")
  out_svt <- PsiNorm(svt)
  expect_true(is(out_svt, "SVT_SparseMatrix"))
  expect_equal(out1, as.matrix(out_svt))

  dm1 <- DelayedArray::DelayedArray(m)
  out_dm1 <- PsiNorm(dm1)
  expect_true(is(out_dm1, "DelayedMatrix"))
  expect_equal(out1, as.matrix(out_dm1))

  dm2 <- DelayedArray::DelayedArray(svt)
  out_dm2 <- PsiNorm(dm2)
  expect_true(is(out_dm2, "DelayedMatrix"))
  expect_equal(out1, as.matrix(out_dm2))

  se <- SummarizedExperiment(m)
  sce <- SingleCellExperiment(m)

  out_se <- PsiNorm(se)
  expect_equal(out1, assay(out_se, "PsiNorm"))

  out_sce <- PsiNorm(sce, whichAssay = 1)
  out2 <- t(t(m)/sizeFactors(out_sce))
  expect_equal(out1, out2)
})
