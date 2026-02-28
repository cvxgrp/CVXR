# Extracted from test-cvxpy-parity.R:176

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
expr <- abs(-x^2)
expect_equal(expr_curvature(expr), CONVEX)
