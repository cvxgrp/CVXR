# Extracted from test-cvxpy-atoms-parity.R:340

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_true(is_nonneg(SumEntries(Constant(1))))
