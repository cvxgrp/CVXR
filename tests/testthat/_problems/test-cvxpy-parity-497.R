# Extracted from test-cvxpy-parity.R:497

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Minimum(Constant(0), Constant(0))), ZERO_SIGN)
