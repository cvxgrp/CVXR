# Extracted from test-cvxpy-atoms-parity.R:1147

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
A <- Variable(c(3L, 3L))
expr <- A[1:2, 1:2]
expect_equal(expr_curvature(expr), AFFINE)
