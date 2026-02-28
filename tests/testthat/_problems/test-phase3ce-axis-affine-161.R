# Extracted from test-phase3ce-axis-affine.R:161

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
n <- NormInf(x)
