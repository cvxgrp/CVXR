# Extracted from test-phase2c-operators.R:505

# test -------------------------------------------------------------------------
d <- DivExpression(Variable(c(3, 1)), Constant(2))
expect_equal(d@shape, c(3L, 1L))
expect_true(S7_inherits(d@args[[2L]], Promote))
