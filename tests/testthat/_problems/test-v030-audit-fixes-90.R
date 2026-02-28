# Extracted from test-v030-audit-fixes.R:90

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
p <- Variable(3, pos = TRUE)
s <- sum_entries(p)
expect_true(is_atom_log_log_convex(s))
