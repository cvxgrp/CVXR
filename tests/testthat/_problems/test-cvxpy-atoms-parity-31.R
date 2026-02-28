# Extracted from test-cvxpy-atoms-parity.R:31

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
y <- Variable(2)
expr <- x + y
atom <- norm_inf(expr)
expect_equal(atom@shape, c(1L, 1L))
expect_equal(expr_curvature(atom), CONVEX)
