# Extracted from test-cvxpy-atoms-parity.R:947

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
e <- Entr(x)
