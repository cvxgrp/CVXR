# Extracted from test-phase3ce-axis-affine.R:280

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 3))
value(x) <- diag(c(1, 2, 3))
tr <- Trace(x)
