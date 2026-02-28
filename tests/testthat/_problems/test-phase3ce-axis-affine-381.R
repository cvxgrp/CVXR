# Extracted from test-phase3ce-axis-affine.R:381

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 3))
u <- upper_tri(x)
expect_true(S7_inherits(u, UpperTri))
