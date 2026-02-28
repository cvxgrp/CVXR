# Extracted from test-phase2d-dispatch.R:168

# test -------------------------------------------------------------------------
x <- Variable(c(3, 4))
A <- matrix(1:8, 4, 2)
result <- x %*% A
expect_true(S7_inherits(result, MulExpression))
