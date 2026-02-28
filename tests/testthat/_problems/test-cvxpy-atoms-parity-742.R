# Extracted from test-cvxpy-atoms-parity.R:742

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(as.numeric(value(Huber(Constant(2), M = 1))), 3, tolerance = 1e-10)
