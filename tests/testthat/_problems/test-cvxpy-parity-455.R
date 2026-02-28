# Extracted from test-cvxpy-parity.R:455

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Maximum(Constant(0), Constant(0))), ZERO_SIGN)
