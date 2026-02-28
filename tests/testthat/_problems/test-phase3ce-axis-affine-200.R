# Extracted from test-phase3ce-axis-affine.R:200

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
n <- cvxr_norm(x, p = 1)
expect_true(S7_inherits(n, Norm1))
