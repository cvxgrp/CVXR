# Extracted from test-phase3ce-axis-affine.R:349

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
A <- matrix(1:4, 2, 2)
x <- Variable(c(3, 3))
k <- Kron(A, x)
