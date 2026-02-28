# Extracted from test-cvxpy-parity.R:447

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Maximum(Variable(1), Variable(1))), UNKNOWN_SIGN)
