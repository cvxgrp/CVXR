# Extracted from test-cvxpy-atoms-parity.R:260

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(MinEntries(Constant(1))), NONNEG_SIGN)
