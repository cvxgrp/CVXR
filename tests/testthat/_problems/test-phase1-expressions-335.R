# Extracted from test-phase1-expressions.R:335

# test -------------------------------------------------------------------------
CVXR:::reset_expr_id()
p <- Parameter(c(2, 3))
expect_equal(p@shape, c(2L, 3L))
expect_equal(expr_name(p), "param1")
