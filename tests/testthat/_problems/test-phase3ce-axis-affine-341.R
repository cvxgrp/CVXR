# Extracted from test-phase3ce-axis-affine.R:341

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
y <- Variable(c(4, 3))
v <- vstack(x, y)
expect_true(S7_inherits(v, VStack))
