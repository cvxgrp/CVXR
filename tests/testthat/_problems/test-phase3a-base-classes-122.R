# Extracted from test-phase3a-base-classes.R:122

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 4))
a <- AxisAtom(x, axis = 1L, keepdims = TRUE)
