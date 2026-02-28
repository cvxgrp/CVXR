# Extracted from test-cvxpy-atoms-parity.R:586

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
A <- Variable(c(2L, 2L))
expr <- DiagMat(A)
expect_equal(expr_sign_str(expr), UNKNOWN_SIGN)
