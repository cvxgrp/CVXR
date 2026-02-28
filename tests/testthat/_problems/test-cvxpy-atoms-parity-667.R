# Extracted from test-cvxpy-atoms-parity.R:667

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
A <- Variable(c(4L, 4L))
ut <- UpperTri(A)
