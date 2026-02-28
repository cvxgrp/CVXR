# Extracted from test-cvxpy-atoms-parity.R:577

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
expr <- DiagVec(x)
expect_equal(expr_sign_str(expr), UNKNOWN_SIGN)
