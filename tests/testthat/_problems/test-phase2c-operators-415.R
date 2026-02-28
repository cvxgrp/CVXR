# Extracted from test-phase2c-operators.R:415

# test -------------------------------------------------------------------------
m <- Multiply(Constant(2), Variable(3))
expect_true(S7_inherits(m, Multiply))
