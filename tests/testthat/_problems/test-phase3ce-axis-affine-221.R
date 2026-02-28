# Extracted from test-phase3ce-axis-affine.R:221

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
r <- Reshape(x, c(6L, 1L))
expect_equal(r@shape, c(6L, 1L))
expect_true(S7_inherits(r, AffAtom))
