# Extracted from test-phase3d-nonaffine.R:43

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(1, nonneg = TRUE)
q <- QuadOverLin(x, y)
