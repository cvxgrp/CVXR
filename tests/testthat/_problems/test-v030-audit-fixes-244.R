# Extracted from test-v030-audit-fixes.R:244

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
dv <- DiagVec(x, k = 2L)
expect_equal(get_data(dv), list(2L))
