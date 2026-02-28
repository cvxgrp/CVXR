# Extracted from test-phase3ce-axis-affine.R:132

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
n <- Norm1(x)
