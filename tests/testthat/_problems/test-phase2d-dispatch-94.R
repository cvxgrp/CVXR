# Extracted from test-phase2d-dispatch.R:94

# test -------------------------------------------------------------------------
x <- Variable(3)
result <- x - 1
expect_true(S7_inherits(result, AddExpression))
