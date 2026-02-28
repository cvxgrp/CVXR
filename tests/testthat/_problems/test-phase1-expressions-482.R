# Extracted from test-phase1-expressions.R:482

# test -------------------------------------------------------------------------
x <- Variable(3)
expect_equal(expr_curvature(x), "AFFINE")
