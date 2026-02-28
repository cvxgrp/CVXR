# Extracted from test-phase3ce-axis-affine.R:36

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
value(x) <- matrix(1:6, 2, 3)
s <- SumEntries(x)
