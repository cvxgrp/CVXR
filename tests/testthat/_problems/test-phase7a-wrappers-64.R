# Extracted from test-phase7a-wrappers.R:64

# test -------------------------------------------------------------------------
x <- Variable(c(5L, 1L))
value(x) <- matrix(c(1, 3, 6, 10, 15), 5, 1)
d <- cvxr_diff(x)
expect_true(S7_inherits(d, Expression))
