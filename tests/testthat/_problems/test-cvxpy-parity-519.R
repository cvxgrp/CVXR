# Extracted from test-cvxpy-parity.R:519

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(MaxEntries(Variable(1))), UNKNOWN_SIGN)
