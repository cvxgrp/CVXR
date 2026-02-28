# Extracted from test-cvxpy-parity.R:527

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(MinEntries(Constant(1))), NONNEG_SIGN)
