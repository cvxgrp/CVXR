# Extracted from test-phase3d-nonaffine.R:35

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
P <- diag(3)
q <- quad_form(x, P)
expect_true(S7_inherits(q, QuadForm))
