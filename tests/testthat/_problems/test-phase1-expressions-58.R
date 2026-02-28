# Extracted from test-phase1-expressions.R:58

# test -------------------------------------------------------------------------
x <- Variable(3)
expect_false(is_nonneg(x))
expect_false(is_nonpos(x))
expect_false(is_zero(x))
expect_equal(expr_sign(x), UNKNOWN_SIGN)
