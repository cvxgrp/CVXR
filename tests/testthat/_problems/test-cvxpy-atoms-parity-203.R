# Extracted from test-cvxpy-atoms-parity.R:203

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
a <- Variable(1)
atom <- quad_over_lin(x^2, a)
expect_equal(expr_curvature(atom), CONVEX)
