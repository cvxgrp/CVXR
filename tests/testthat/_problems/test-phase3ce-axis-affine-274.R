# Extracted from test-phase3ce-axis-affine.R:274

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 3))
tr <- Trace(x)
expect_equal(tr@shape, c(1L, 1L))
expect_true(S7_inherits(tr, AffAtom))
