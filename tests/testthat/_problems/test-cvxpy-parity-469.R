# Extracted from test-cvxpy-parity.R:469

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Minimum(Constant(1), Constant(2))), NONNEG_SIGN)
