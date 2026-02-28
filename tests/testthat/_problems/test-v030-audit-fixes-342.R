# Extracted from test-v030-audit-fixes.R:342

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
q <- quad_over_lin(x, Constant(1))
expect_true(is_atom_log_log_convex(q))
