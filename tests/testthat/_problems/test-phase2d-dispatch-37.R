# Extracted from test-phase2d-dispatch.R:37

# test -------------------------------------------------------------------------
x <- Variable(3)
A <- matrix(1:3, 3, 1)
result <- A + x
expect_true(S7_inherits(result, AddExpression))
