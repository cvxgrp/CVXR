# Extracted from test-phase3ce-axis-affine.R:82

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 4))
m <- MaxEntries(x, axis = 2L)
