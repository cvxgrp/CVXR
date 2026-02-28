# Extracted from test-cvxpy-parity.R:259

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
y <- Variable(3)
A <- Constant(matrix(rnorm(12), 4, 3))
b <- Constant(matrix(rnorm(4), ncol = 1))
Aeq <- Constant(matrix(rnorm(6), 2, 3))
beq <- Constant(matrix(rnorm(2), ncol = 1))
p <- Problem(Minimize(sum_squares(A %*% y - b)),
               list(Maximum(Constant(matrix(1, 3, 1)), 3 * y) <= 200,
                    Abs(2 * y) <= 100,
                    norm1(2 * y) <= 1000,
                    Aeq %*% y == beq))
