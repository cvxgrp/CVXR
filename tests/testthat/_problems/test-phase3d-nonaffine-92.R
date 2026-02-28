# Extracted from test-phase3d-nonaffine.R:92

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(5)
value(x) <- c(1, 5, 2, 4, 3)
s <- SumLargest(x, k = 2L)
