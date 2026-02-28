# Extracted from test-phase2d-dispatch.R:297

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(3)
result <- x - 2 * y
expect_true(S7_inherits(result, AddExpression))
