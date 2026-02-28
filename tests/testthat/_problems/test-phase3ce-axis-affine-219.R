# Extracted from test-phase3ce-axis-affine.R:219

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
r <- Reshape(x, c(6L, 1L))
