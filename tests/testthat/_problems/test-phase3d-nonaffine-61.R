# Extracted from test-phase3d-nonaffine.R:61

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(1)
q <- QuadOverLin(x, y)
