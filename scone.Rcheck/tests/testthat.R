Sys.setenv("R_TESTS" = "")

library(testthat)
library(scone)

test_check("scone")
