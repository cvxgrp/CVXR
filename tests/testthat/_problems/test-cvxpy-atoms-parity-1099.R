# Extracted from test-cvxpy-atoms-parity.R:1099

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
a <- Constant(matrix(c(1, 1, 1), 3, 1))
b <- Variable(c(2L, 1L))
value(b) <- matrix(c(1, 2), 2, 1)
expr <- Convolve(a, b)
