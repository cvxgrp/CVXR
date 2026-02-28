# Extracted from test-cvxpy-atoms-parity.R:302

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expr <- Maximum(Constant(1), Variable(2))
