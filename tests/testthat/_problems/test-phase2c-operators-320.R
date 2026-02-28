# Extracted from test-phase2c-operators.R:320

# test -------------------------------------------------------------------------
m <- MulExpression(Constant(matrix(1:6, 3, 2)), Variable(c(2, 1)))
expect_true(S7_inherits(m, MulExpression))
