# Extracted from test-phase1-expressions.R:204

# test -------------------------------------------------------------------------
c1 <- Constant(5)
expect_true(is_constant(c1))
expect_true(is_affine(c1))
expect_true(is_convex(c1))
expect_true(is_concave(c1))
expect_true(is_dcp(c1))
expect_equal(expr_curvature(c1), "CONSTANT")
