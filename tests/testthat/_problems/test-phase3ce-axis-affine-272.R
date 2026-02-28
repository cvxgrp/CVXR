# Extracted from test-phase3ce-axis-affine.R:272

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 3))
tr <- Trace(x)
