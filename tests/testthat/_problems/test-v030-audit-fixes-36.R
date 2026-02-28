# Extracted from test-v030-audit-fixes.R:36

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_true(is_pos(Constant(5)))
