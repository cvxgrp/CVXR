# Extracted from test-cvxpy-atoms-parity.R:507

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
A <- Variable(c(2L, 2L))
expr <- reshape_expr(A, c(4L, 1L), order = "F")
expect_equal(expr_sign_str(expr), UNKNOWN_SIGN)
