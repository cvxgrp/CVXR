# Extracted from test-cvxpy-parity.R:170

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
expr <- abs(1 + exp(x))
expect_equal(expr_curvature(expr), CONVEX)
