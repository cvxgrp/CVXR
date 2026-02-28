# Extracted from test-phase1-expressions.R:37

# test -------------------------------------------------------------------------
x <- Variable(3, var_id = 999L)
expect_equal(x@id, 999L)
expect_equal(expr_name(x), "var999")
