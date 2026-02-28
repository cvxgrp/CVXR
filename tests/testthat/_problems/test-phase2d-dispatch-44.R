# Extracted from test-phase2d-dispatch.R:44

# test -------------------------------------------------------------------------
a <- Constant(5)
b <- Constant(3)
result <- a + b
expect_true(S7_inherits(result, AddExpression))
