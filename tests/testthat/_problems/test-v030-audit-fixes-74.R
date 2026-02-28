# Extracted from test-v030-audit-fixes.R:74

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
p <- Variable(3, pos = TRUE)
e <- exp(p)
expect_true(is_atom_log_log_convex(e))
