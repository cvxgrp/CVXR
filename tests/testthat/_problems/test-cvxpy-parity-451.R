# Extracted from test-cvxpy-parity.R:451

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Maximum(Variable(1), Constant(-2))), UNKNOWN_SIGN)
