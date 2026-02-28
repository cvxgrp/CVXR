# Extracted from test-v030-audit-fixes.R:53

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
p <- Variable(3, pos = TRUE)
expect_true(is_log_log_convex(p))
