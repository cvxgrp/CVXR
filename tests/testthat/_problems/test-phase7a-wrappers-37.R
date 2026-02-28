# Extracted from test-phase7a-wrappers.R:37

# test -------------------------------------------------------------------------
x <- Variable(c(5L, 1L))
d <- cvxr_diff(x, k = 0L)
expect_true(S7_inherits(d, Expression))
