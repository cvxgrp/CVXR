# Extracted from test-cvxpy-parity.R:489

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Minimum(Variable(1), Variable(1))), UNKNOWN_SIGN)
