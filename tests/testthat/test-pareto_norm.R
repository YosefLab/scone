test_that("PSiNorm works with all input classes", {
  m <- matrix(c(1,0,2,0,2,9,3,0), ncol=2)
  
  out1 <- PsiNorm(m)
  
  se <- SummarizedExperiment(m)
  sce <- SingleCellExperiment(m)
  
  out_se <- PsiNorm(se)
  expect_equal(out1, assay(out_se, "PsiNorm"))
  
  out_sce <- PsiNorm(sce, whichAssay = 1)
  out2 <- t(t(m)/sizeFactors(out_sce))
  expect_equal(out1, out2)
})
