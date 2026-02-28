# Extracted from test-cvxpy-parity.R:120

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
expr <- 1 + exp(x)
expect_equal(expr_curvature(expr), CONVEX)
