# Extracted from test-v030-audit-fixes.R:168

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
q <- quad_over_lin(x, Constant(1), axis = 2L, keepdims = TRUE)
d <- get_data(q)
