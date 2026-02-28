# Extracted from test-cvxpy-atoms-parity.R:486

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
A <- Variable(c(2L, 2L))
B <- Variable(c(2L, 2L))
atom <- HStack(A, B)
