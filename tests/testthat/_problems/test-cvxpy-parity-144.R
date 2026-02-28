# Extracted from test-cvxpy-parity.R:144

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
expr <- -(exp(x))^2
expect_equal(expr_curvature(expr), CONCAVE)
