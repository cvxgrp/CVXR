# Extracted from test-cvxpy-atoms-parity.R:662

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
C <- Variable(c(3L, 2L))
expect_error(UpperTri(C), "square")
