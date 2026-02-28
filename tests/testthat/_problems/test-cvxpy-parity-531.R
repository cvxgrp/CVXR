# Extracted from test-cvxpy-parity.R:531

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(MinEntries(Constant(-2))), NONPOS_SIGN)
