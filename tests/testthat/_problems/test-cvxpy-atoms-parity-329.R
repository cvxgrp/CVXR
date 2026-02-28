# Extracted from test-cvxpy-atoms-parity.R:329

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expr <- Minimum(Constant(-1), Variable(2))
