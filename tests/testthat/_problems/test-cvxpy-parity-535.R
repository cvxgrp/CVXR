# Extracted from test-cvxpy-parity.R:535

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(MinEntries(Variable(1))), UNKNOWN_SIGN)
