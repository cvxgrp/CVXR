# Extracted from test-cvxpy-atoms-parity.R:521

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
expr <- reshape_expr(x^2, c(1L, 2L), order = "F")
expect_true(is_nonneg(expr))
expect_equal(expr_curvature(expr), CONVEX)
