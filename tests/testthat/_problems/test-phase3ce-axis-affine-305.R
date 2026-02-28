# Extracted from test-phase3ce-axis-affine.R:305

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 2))
y <- Variable(c(2, 3))
value(x) <- matrix(1:4, 2, 2)
value(y) <- matrix(5:10, 2, 3)
h <- HStack(x, y)
