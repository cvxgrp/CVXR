# Extracted from test-phase3a-base-classes.R:57

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(3)
e <- Elementwise(args = list(x, y))
