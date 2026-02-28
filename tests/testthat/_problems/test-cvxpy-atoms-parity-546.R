# Extracted from test-cvxpy-atoms-parity.R:546

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
C <- Variable(c(3L, 2L))
expr <- vec(C)
expect_equal(expr_sign_str(expr), UNKNOWN_SIGN)
