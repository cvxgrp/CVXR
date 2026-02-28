# Extracted from test-phase2d-dispatch.R:290

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(3)
result <- (x + y) / 2
expect_true(S7_inherits(result, DivExpression))
