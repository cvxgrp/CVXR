# Extracted from test-cvxpy-parity.R:473

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Minimum(Constant(1), Variable(1))), UNKNOWN_SIGN)
