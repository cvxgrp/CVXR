# Extracted from test-phase1-expressions.R:10

# test -------------------------------------------------------------------------
CVXR:::reset_expr_id()
x <- Variable(3)
expect_equal(x@shape, c(3L, 1L))
expect_equal(expr_name(x), "var1")
