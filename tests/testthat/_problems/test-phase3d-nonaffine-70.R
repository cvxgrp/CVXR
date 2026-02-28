# Extracted from test-phase3d-nonaffine.R:70

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(1, nonneg = TRUE)
q <- quad_over_lin(x, y)
expect_true(S7_inherits(q, QuadOverLin))
