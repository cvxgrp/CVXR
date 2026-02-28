# Extracted from test-phase2c-operators.R:86

# test -------------------------------------------------------------------------
x <- Constant(5)
y <- Variable(c(3, 1))
result <- broadcast_args(x, y)
expect_true(S7_inherits(result[[1L]], Promote))
