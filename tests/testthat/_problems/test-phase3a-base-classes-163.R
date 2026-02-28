# Extracted from test-phase3a-base-classes.R:163

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 4))
a <- AxisAffAtom(x, axis = 2L)
