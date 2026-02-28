# Extracted from test-cvxpy-atoms-parity.R:451

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
y <- Variable(2)
A <- Variable(c(2L, 2L))
C <- Variable(c(3L, 2L))
B <- Variable(c(2L, 2L))
atom <- VStack(x, y, x)
