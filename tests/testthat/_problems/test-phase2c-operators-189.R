# Extracted from test-phase2c-operators.R:189

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(3)
a <- AddExpression(list(x, y))
expect_true(S7_inherits(a, AddExpression))
