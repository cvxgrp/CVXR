# Extracted from test-phase2c-operators.R:113

# test -------------------------------------------------------------------------
n <- NegExpression(Variable(3))
expect_true(S7_inherits(n, NegExpression))
