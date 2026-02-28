# Extracted from test-phase3a-base-classes.R:76

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
e <- Elementwise(args = list(x))
