# Extracted from test-cvxpy-parity.R:463

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Maximum(Constant(-3), Constant(-2))), NONPOS_SIGN)
