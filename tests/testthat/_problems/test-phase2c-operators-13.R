# Extracted from test-phase2c-operators.R:13

# test -------------------------------------------------------------------------
expect_true(inherits(Promote, "S7_class"))
p <- Promote(Constant(1), c(3L, 1L))
expect_true(S7_inherits(p, Promote))
