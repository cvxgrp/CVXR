# Extracted from test-phase3ce-axis-affine.R:332

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 2))
y <- Variable(c(3, 2))
value(x) <- matrix(1:4, 2, 2)
value(y) <- matrix(5:10, 3, 2)
v <- VStack(x, y)
