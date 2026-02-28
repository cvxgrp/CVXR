# Extracted from test-phase3d-nonaffine.R:19

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
P <- diag(3)
value(x) <- c(1, 2, 3)
q <- QuadForm(x, P)
