# Extracted from test-phase3a-base-classes.R:77

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
e <- Elementwise(args = list(x))
expect_true(S7_inherits(e, Atom))
