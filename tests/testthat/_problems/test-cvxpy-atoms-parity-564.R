# Extracted from test-cvxpy-atoms-parity.R:564

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
a <- Variable(1)
expr <- vec(a^2)
expect_true(is_nonneg(expr))
expect_equal(expr_curvature(expr), CONVEX)
