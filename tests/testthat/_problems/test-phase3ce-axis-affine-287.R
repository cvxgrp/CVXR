# Extracted from test-phase3ce-axis-affine.R:287

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(3, 3))
tr <- matrix_trace(x)
expect_true(S7_inherits(tr, Trace))
