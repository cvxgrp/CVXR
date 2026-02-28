# Extracted from test-cvxpy-atoms-parity.R:757

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
atom <- SumLargest(x, k = 2)
