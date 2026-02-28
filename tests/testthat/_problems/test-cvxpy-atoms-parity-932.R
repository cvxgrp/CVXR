# Extracted from test-cvxpy-atoms-parity.R:932

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
e <- Exp(x)
