# Extracted from test-phase2d-dispatch.R:429

# test -------------------------------------------------------------------------
x <- Variable(c(3, 4))
result <- 2 * x
expect_true(S7_inherits(result, Multiply))
