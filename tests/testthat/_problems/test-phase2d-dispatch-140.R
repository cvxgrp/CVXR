# Extracted from test-phase2d-dispatch.R:140

# test -------------------------------------------------------------------------
x <- Variable(3)
result <- x / 2
expect_true(S7_inherits(result, DivExpression))
