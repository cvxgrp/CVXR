# Extracted from test-cvxpy-parity.R:431

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Maximum(Constant(1), Variable(1))), NONNEG_SIGN)
