# Extracted from test-cvxpy-parity.R:477

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Minimum(Constant(1), Constant(-2))), NONPOS_SIGN)
