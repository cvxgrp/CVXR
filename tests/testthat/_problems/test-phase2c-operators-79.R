# Extracted from test-phase2c-operators.R:79

# test -------------------------------------------------------------------------
x <- Variable(c(3, 1))
expect_error(cvxr_promote(x, c(4L, 1L)), "Only scalars")
