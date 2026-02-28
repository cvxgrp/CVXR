# Extracted from test-cvxpy-atoms-parity.R:787

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
atom <- sum_smallest(x, 2)
expect_true(is_pwl(atom))
