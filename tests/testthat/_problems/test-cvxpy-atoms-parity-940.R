# Extracted from test-cvxpy-atoms-parity.R:940

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
l <- Log(x)
