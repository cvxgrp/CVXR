# Extracted from test-phase3ce-axis-affine.R:297

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 2))
y <- Variable(c(3, 4))
h <- HStack(x, y)
expect_equal(h@shape, c(3L, 6L))
expect_true(S7_inherits(h, AffAtom))
