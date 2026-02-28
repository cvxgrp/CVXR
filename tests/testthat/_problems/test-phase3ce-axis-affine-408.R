# Extracted from test-phase3ce-axis-affine.R:408

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(5)
c_atom <- Cumsum(x)
expect_equal(c_atom@shape, c(5L, 1L))
expect_true(S7_inherits(c_atom, AxisAffAtom))
