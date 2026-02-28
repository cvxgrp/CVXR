# Extracted from test-v030-audit-fixes.R:13

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
c1 <- Constant(5)
expect_true(is_log_log_constant(c1))
