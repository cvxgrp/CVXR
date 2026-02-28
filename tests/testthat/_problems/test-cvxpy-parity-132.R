# Extracted from test-cvxpy-parity.R:132

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
expr <- 1 - sqrt(x)
expect_equal(expr_curvature(expr), CONVEX)
