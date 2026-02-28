# Extracted from test-phase2d-dispatch.R:406

# test -------------------------------------------------------------------------
p <- Parameter(3)
x <- Variable(3)
result <- p * x
expect_true(S7_inherits(result, Multiply))
