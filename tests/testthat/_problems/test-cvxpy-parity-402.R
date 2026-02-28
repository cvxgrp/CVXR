# Extracted from test-cvxpy-parity.R:402

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expr_nonneg <- Variable(1, nonneg = TRUE)
tmp_nn <- Minimum(expr_nonneg, Constant(Inf))
