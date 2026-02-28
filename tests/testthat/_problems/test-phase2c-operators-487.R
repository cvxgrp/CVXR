# Extracted from test-phase2c-operators.R:487

# test -------------------------------------------------------------------------
d <- DivExpression(Variable(3), Constant(2))
expect_true(S7_inherits(d, DivExpression))
