# Extracted from test-phase3ce-axis-affine.R:154

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
n <- norm1(x)
expect_true(S7_inherits(n, Norm1))
