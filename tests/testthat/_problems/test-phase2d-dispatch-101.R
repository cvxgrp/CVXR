# Extracted from test-phase2d-dispatch.R:101

# test -------------------------------------------------------------------------
x <- Variable(3)
result <- 1 - x
expect_true(S7_inherits(result, AddExpression))
