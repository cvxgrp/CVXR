# Extracted from test-cvxpy-parity.R:511

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(MaxEntries(Constant(1))), NONNEG_SIGN)
