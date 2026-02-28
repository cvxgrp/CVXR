# Extracted from test-cvxpy-atoms-parity.R:1089

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
a <- Constant(matrix(c(1, 0, 0, 1), 2, 2))
b <- Constant(matrix(c(1, 2, 3, 4), 2, 2))
k <- Kron(a, b)
