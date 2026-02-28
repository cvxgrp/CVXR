# Extracted from test-cvxpy-atoms-parity.R:234

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(max_entries(Constant(1))), NONNEG_SIGN)
