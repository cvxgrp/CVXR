# Extracted from test-v030-audit-fixes.R:151

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(c(2, 3))
xval <- matrix(1:6, nrow = 2, ncol = 3)
q_null <- quad_over_lin(x, Constant(2))
q_ax0 <- quad_over_lin(x, Constant(2), axis = 2L)
q_ax1 <- quad_over_lin(x, Constant(2), axis = 1L)
nv_null <- numeric_value(q_null, list(xval, matrix(2, 1, 1)))
