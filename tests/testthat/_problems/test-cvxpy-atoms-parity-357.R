# Extracted from test-cvxpy-atoms-parity.R:357

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(SumEntries(Variable(c(2L, 1L)), keepdims = TRUE)@shape, c(1L, 1L))
