# Extracted from test-phase2d-dispatch.R:128

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(3)
result <- x * y
expect_true(S7_inherits(result, Multiply))
