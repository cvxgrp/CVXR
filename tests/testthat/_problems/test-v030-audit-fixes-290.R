# Extracted from test-v030-audit-fixes.R:290

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(4, 4))
dm <- DiagMat(x, k = -1L)
expect_equal(get_data(dm), list(-1L))
