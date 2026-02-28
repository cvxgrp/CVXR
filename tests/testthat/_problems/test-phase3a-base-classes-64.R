# Extracted from test-phase3a-base-classes.R:64

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
c1 <- Constant(1)
e <- Elementwise(args = list(x, c1))
