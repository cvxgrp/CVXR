# Extracted from test-cvxpy-atoms-parity.R:824

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
a <- Constant(matrix(1, 3, 1))
b_nonneg <- Constant(matrix(c(1, 2), 2, 1))
expr <- Convolve(a, b_nonneg)
