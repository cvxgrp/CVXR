# Extracted from test-cvxpy-atoms-parity.R:362

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(SumEntries(Variable(2), axis = 2L)@shape, c(1L, 1L))
