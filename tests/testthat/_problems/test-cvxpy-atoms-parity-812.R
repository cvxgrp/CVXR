# Extracted from test-cvxpy-atoms-parity.R:812

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
expect_error(Kron(x, x), "constant")
