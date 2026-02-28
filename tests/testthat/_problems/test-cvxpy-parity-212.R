# Extracted from test-cvxpy-parity.R:212

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
expr <- abs(entr(x))
expect_equal(expr_curvature(expr), UNKNOWN_CURVATURE)
