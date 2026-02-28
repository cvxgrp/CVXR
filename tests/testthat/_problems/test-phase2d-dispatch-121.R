# Extracted from test-phase2d-dispatch.R:121

# test -------------------------------------------------------------------------
x <- Variable(3)
result <- 2 * x
expect_true(S7_inherits(result, Multiply))
