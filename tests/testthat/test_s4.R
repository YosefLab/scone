context("Tests for the SconeExperiment object")
set.seed(13124)

nrows <- 200
ncols <- 6
counts <- matrix(rpois(nrows * ncols, lambda=10), nrows)
rowdata <- data.frame(poscon=c(rep(TRUE, 10), rep(FALSE, nrows-10)))
coldata <- data.frame(bio=gl(2, 3))
se <- SummarizedExperiment(assays=SimpleList(counts=counts),
                           rowData=rowdata, colData=coldata)

test_that("The two constructors are equivalent", {
  expect_equal(sconeExperiment(assay(se)), sconeExperiment(assay(se)))

  scone1 <- sconeExperiment(assay(se), bio=coldata$bio, poscon=rowdata$poscon)
  scone2 <- sconeExperiment(se, which_bio=1L, which_poscon=1L)

  expect_equal(scone1, scone2)
  }
)
