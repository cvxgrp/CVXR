# Extracted from test-cvxpy-atoms-parity.R:924

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
a <- Abs(x)
