# Extracted from test-phase3ce-axis-affine.R:12

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 4))
s <- SumEntries(x)
expect_equal(s@shape, c(1L, 1L))
expect_true(S7_inherits(s, AxisAffAtom))
