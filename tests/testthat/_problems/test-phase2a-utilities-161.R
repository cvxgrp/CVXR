# Extracted from test-phase2a-utilities.R:161

# test -------------------------------------------------------------------------
expect_error(sum_shapes(list(c(3L, 4L), c(5L, 4L))),
               "Cannot broadcast")
