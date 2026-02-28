# Extracted from test-v030-audit-fixes.R:109

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
p <- Variable(3, pos = TRUE)
m <- min_entries(p)
expect_false(is_atom_log_log_convex(m))
