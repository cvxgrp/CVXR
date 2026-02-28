# Extracted from test-phase2d-dispatch.R:175

# test -------------------------------------------------------------------------
x <- Variable(c(3, 4))
y <- Variable(c(4, 2))
result <- x %*% y
expect_true(S7_inherits(result, MulExpression))
