# Extracted from test-phase3ce-axis-affine.R:177

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
n <- Pnorm(x, p = 2)
