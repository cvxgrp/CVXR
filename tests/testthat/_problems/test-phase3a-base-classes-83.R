# Extracted from test-phase3a-base-classes.R:83

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 3))
e <- Elementwise(args = list(x))
