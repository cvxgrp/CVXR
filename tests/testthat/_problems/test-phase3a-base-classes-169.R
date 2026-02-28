# Extracted from test-phase3a-base-classes.R:169

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
a <- AxisAffAtom(x)
