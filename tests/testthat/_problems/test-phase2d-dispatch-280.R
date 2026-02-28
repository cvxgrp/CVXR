# Extracted from test-phase2d-dispatch.R:280

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(3)
result <- 2 * x + 3 * y
expect_true(S7_inherits(result, AddExpression))
