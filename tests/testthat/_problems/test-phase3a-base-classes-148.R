# Extracted from test-phase3a-base-classes.R:148

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
a <- AxisAtom(x)
