# Extracted from test-cvxpy-atoms-parity.R:81

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2L, 1L))
y <- Variable(c(2L, 1L))
expr <- x + y
expect_equal(expr_curvature(power(expr, 2)), CONVEX)
