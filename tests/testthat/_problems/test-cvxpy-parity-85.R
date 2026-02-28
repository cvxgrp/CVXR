# Extracted from test-cvxpy-parity.R:85

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
cvx <- Variable(1)^2
aff <- Variable(1)
expect_equal(expr_curvature(-cvx), CONCAVE)
