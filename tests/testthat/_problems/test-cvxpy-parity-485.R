# Extracted from test-cvxpy-parity.R:485

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Minimum(Variable(1), Constant(0))), NONPOS_SIGN)
