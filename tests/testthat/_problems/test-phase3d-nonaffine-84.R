# Extracted from test-phase3d-nonaffine.R:84

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(5)
s <- SumLargest(x, k = 3L)
