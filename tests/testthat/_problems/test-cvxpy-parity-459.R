# Extracted from test-cvxpy-parity.R:459

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Maximum(Constant(0), Constant(-2))), ZERO_SIGN)
