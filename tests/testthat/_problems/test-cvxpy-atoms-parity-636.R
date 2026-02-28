# Extracted from test-cvxpy-atoms-parity.R:636

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
A <- Variable(c(2L, 2L))
expr <- Trace(A)
