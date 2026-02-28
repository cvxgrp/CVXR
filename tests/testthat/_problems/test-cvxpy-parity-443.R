# Extracted from test-cvxpy-parity.R:443

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
expect_equal(expr_sign_str(Maximum(Variable(1), Constant(0))), NONNEG_SIGN)
