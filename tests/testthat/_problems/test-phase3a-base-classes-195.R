# Extracted from test-phase3a-base-classes.R:195

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
a <- Elementwise(args = list(x))
