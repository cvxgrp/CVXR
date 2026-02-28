# Extracted from test-cvxpy-atoms-parity.R:1039

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
v <- c(10, 1, 5, 3, 7)
x <- Constant(matrix(v, ncol = 1))
expect_equal(as.numeric(value(SumLargest(x, k = 2))), 17, tolerance = 1e-10)
