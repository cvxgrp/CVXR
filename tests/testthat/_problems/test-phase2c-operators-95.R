# Extracted from test-phase2c-operators.R:95

# test -------------------------------------------------------------------------
x <- Variable(c(3, 1))
y <- Constant(5)
result <- broadcast_args(x, y)
expect_identical(result[[1L]], x)
expect_true(S7_inherits(result[[2L]], Promote))
