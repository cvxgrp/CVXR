# Extracted from test-cvxpy-parity.R:269

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
y <- Variable(3)
A <- Constant(matrix(rnorm(12), 4, 3))
b <- Constant(matrix(rnorm(4), ncol = 1))
p <- Problem(Minimize(sum_squares(A %*% y - b)),
               list(Maximum(Constant(matrix(1, 3, 1)), 3 * y^2) <= 200))
