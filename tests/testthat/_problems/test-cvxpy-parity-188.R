# Extracted from test-cvxpy-parity.R:188

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1, nonneg = TRUE)
expr <- abs(-x)
expect_equal(expr_curvature(expr), CONVEX)
