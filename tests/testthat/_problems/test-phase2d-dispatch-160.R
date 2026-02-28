# Extracted from test-phase2d-dispatch.R:160

# test -------------------------------------------------------------------------
A <- matrix(1:6, 3, 2)
x <- Variable(c(2, 1))
result <- A %*% x
expect_true(S7_inherits(result, MulExpression))
