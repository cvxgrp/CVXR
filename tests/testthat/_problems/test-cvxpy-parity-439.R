# Extracted from test-cvxpy-parity.R:439

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Maximum(Constant(1), Constant(0))), NONNEG_SIGN)
