# Extracted from test-cvxpy-parity.R:380

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expr_nonneg <- Variable(1, nonneg = TRUE)
tmp_nn <- Maximum(expr_nonneg, Constant(-Inf))
