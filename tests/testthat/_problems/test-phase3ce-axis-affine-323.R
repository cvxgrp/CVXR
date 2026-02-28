# Extracted from test-phase3ce-axis-affine.R:323

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
y <- Variable(c(4, 3))
v <- VStack(x, y)
