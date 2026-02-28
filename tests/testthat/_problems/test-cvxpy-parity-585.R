# Extracted from test-cvxpy-parity.R:585

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
p <- Problem(Minimize(exp(x)))
d2c <- Dcp2Cone()
