# Extracted from test-cvxpy-atoms-parity.R:839

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
y <- Variable(2)
expect_error(Convolve(x, y), "constant")
