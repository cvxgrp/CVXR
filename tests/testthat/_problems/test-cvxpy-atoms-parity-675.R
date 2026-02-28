# Extracted from test-cvxpy-atoms-parity.R:675

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
A <- Variable(c(3L, 3L))
value(A) <- matrix(1:9, 3, 3)
ut <- UpperTri(A)
