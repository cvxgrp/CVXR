# Extracted from test-phase3ce-axis-affine.R:244

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
d <- DiagVec(x)
expect_equal(d@shape, c(3L, 3L))
expect_true(S7_inherits(d, AffAtom))
