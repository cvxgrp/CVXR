# Extracted from test-phase3ce-axis-affine.R:295

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 2))
y <- Variable(c(3, 4))
h <- HStack(x, y)
