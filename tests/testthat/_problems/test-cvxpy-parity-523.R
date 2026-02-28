# Extracted from test-cvxpy-parity.R:523

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(MaxEntries(Constant(0))), ZERO_SIGN)
