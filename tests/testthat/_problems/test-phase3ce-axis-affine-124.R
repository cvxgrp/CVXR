# Extracted from test-phase3ce-axis-affine.R:124

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
m <- -MinEntries(x)
