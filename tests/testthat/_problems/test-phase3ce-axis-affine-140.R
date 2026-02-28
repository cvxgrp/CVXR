# Extracted from test-phase3ce-axis-affine.R:140

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 4))
n <- Norm1(x, axis = 2L)
