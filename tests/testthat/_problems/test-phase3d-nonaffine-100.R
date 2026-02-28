# Extracted from test-phase3d-nonaffine.R:100

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(5)
s <- sum_largest(x, k = 3L)
expect_true(S7_inherits(s, SumLargest))
