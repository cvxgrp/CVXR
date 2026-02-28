# Extracted from test-cvxpy-atoms-parity.R:125

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
atom <- p_norm(x, p = 1.5)
expect_equal(atom@shape, c(1L, 1L))
expect_equal(expr_curvature(atom), CONVEX)
