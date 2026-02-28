# Extracted from test-v030-audit-fixes.R:81

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
p <- Variable(3, pos = TRUE)
l <- log(p)
expect_false(is_atom_log_log_convex(l))
