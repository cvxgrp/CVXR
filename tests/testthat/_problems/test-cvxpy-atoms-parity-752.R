# Extracted from test-cvxpy-atoms-parity.R:752

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
expect_error(SumLargest(x, k = -1), "positive|[Ss]econd")
