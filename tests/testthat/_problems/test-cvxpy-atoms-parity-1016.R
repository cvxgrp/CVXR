# Extracted from test-cvxpy-atoms-parity.R:1016

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
A <- Variable(c(3L, 3L))
value(A) <- diag(c(1, 2, 3))
expect_equal(as.numeric(value(Trace(A))), 6, tolerance = 1e-10)
