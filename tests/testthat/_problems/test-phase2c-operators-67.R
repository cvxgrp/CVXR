# Extracted from test-phase2c-operators.R:67

# test -------------------------------------------------------------------------
x <- Variable(1)
p <- cvxr_promote(x, c(3L, 1L))
expect_true(S7_inherits(p, Promote))
