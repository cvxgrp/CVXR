# Extracted from test-phase2c-operators.R:432

# test -------------------------------------------------------------------------
lhs <- Constant(2)
rhs <- Variable(c(3, 1))
m <- Multiply(lhs, rhs)
expect_equal(m@shape, c(3L, 1L))
expect_true(S7_inherits(m@args[[1L]], Promote))
