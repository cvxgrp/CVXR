# Extracted from test-phase3ce-axis-affine.R:398

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
a <- c(1, 1, 1)
x <- Variable(5)
c_expr <- conv(a, x)
expect_true(S7_inherits(c_expr, Convolve))
