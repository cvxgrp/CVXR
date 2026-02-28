# Extracted from test-cvxpy-parity.R:505

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Minimum(Constant(-3), Constant(-2))), NONPOS_SIGN)
