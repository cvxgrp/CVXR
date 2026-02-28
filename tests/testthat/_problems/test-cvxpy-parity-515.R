# Extracted from test-cvxpy-parity.R:515

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(MaxEntries(Constant(-2))), NONPOS_SIGN)
