# Extracted from test-phase3ce-axis-affine.R:169

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
value(x) <- c(-1, 2, -3)
n <- NormInf(x)
