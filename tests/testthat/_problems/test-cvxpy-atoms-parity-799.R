# Extracted from test-cvxpy-atoms-parity.R:799

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
a <- Constant(matrix(1, 3, 2))
b_nonneg <- Constant(matrix(c(1, 2), 2, 1))
expr <- Kron(a, b_nonneg)
