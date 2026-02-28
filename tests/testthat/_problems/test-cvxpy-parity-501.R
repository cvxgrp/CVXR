# Extracted from test-cvxpy-parity.R:501

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Minimum(Constant(0), Constant(-2))), NONPOS_SIGN)
