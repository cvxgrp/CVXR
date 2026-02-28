# Extracted from test-phase2d-dispatch.R:55

# test -------------------------------------------------------------------------
x <- Variable(3)
result <- -x
expect_true(S7_inherits(result, NegExpression))
