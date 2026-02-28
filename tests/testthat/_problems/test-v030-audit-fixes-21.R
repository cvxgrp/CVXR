# Extracted from test-v030-audit-fixes.R:21

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
c0 <- Constant(0)
expect_false(is_log_log_constant(c0))
