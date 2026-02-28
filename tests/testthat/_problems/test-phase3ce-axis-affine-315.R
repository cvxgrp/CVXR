# Extracted from test-phase3ce-axis-affine.R:315

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 2))
y <- Variable(c(3, 4))
h <- hstack(x, y)
expect_true(S7_inherits(h, HStack))
