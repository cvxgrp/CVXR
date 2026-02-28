# Extracted from test-cvxpy-parity.R:493

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Minimum(Variable(1), Constant(-2))), NONPOS_SIGN)
