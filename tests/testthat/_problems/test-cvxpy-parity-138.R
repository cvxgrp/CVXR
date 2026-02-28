# Extracted from test-cvxpy-parity.R:138

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
expr <- log(sqrt(x))
expect_equal(expr_curvature(expr), CONCAVE)
