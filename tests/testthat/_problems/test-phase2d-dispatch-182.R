# Extracted from test-phase2d-dispatch.R:182

# test -------------------------------------------------------------------------
x <- Variable(3)
result <- 5 %*% x
expect_true(S7_inherits(result, MulExpression))
