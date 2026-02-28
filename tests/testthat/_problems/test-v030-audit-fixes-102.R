# Extracted from test-v030-audit-fixes.R:102

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
p <- Variable(3, pos = TRUE)
m <- max_entries(p)
expect_true(is_atom_log_log_convex(m))
