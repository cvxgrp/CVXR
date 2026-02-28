# Extracted from test-cvxpy-atoms-parity.R:429

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
neg <- Constant(matrix(c(-1, -2), ncol = 1))
expr <- neg * (x^2)
expect_equal(expr_curvature(expr), CONCAVE)
