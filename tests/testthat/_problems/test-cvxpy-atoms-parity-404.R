# Extracted from test-cvxpy-atoms-parity.R:404

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
expr <- Constant(matrix(c(1, -1), ncol = 1)) * x
expect_equal(expr_sign_str(expr), UNKNOWN_SIGN)
