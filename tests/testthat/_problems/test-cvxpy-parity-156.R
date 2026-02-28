# Extracted from test-cvxpy-parity.R:156

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1, nonneg = TRUE)
expr <- entr(x)
expect_equal(expr_curvature(expr), CONCAVE)
