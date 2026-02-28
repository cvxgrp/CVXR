# Extracted from test-phase3ce-axis-affine.R:192

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
n <- p_norm(x, p = 2)
expect_true(S7_inherits(n, Pnorm))
