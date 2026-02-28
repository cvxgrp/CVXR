# Extracted from test-phase3ce-axis-affine.R:235

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
r <- reshape_expr(x, c(6L, 1L))
expect_true(S7_inherits(r, Reshape))
