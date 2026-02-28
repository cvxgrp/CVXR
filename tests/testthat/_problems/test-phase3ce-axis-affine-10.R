# Extracted from test-phase3ce-axis-affine.R:10

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 4))
s <- SumEntries(x)
