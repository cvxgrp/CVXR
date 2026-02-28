# Extracted from test-cvxpy-parity.R:539

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(MinEntries(Constant(0))), ZERO_SIGN)
