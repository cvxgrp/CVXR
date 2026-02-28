# Extracted from test-cvxpy-atoms-parity.R:597

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
A <- Variable(c(2L, 2L))
expr <- DiagMat(Transpose(A))
