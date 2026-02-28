# Extracted from test-phase3a-base-classes.R:149

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
a <- AxisAtom(x)
expect_true(S7_inherits(a, Atom))
