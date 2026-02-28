# Extracted from test-cvxpy-atoms-parity.R:728

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
h <- Huber(x, M = 1)
