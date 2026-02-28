# Extracted from test-v030-audit-fixes.R:175

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
q <- quad_over_lin(x, Constant(1))
expect_true(has_quadratic_term(q))
