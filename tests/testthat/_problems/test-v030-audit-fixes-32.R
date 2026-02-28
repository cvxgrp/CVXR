# Extracted from test-v030-audit-fixes.R:32

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
expect_false(is_log_log_constant(x))
