# Extracted from test-cvxpy-atoms-parity.R:1060

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(as.numeric(value(Huber(Constant(0.5), M = 1))), 0.25, tolerance = 1e-10)
