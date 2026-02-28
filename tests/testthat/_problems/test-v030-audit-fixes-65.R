# Extracted from test-v030-audit-fixes.R:65

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
c1 <- Constant(3)
expect_true(is_log_log_affine(c1))
