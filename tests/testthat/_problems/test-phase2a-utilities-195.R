# Extracted from test-phase2a-utilities.R:195

# test -------------------------------------------------------------------------
expect_error(mul_shapes(c(3L, 4L), c(5L, 2L)),
               "Incompatible dimensions")
