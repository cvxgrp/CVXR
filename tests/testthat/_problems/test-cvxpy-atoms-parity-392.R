# Extracted from test-cvxpy-atoms-parity.R:392

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
A <- Matrix::Diagonal(3)
expect_equal(as.numeric(value(SumEntries(Constant(A)))), 3)
