# Extracted from test-cvxpy-parity.R:194

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
expr <- abs(x)
expect_equal(expr_curvature(expr), CONVEX)
