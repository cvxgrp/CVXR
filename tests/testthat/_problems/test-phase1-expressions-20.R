# Extracted from test-phase1-expressions.R:20

# test -------------------------------------------------------------------------
x <- Variable(3, name = "myvar")
expect_equal(expr_name(x), "myvar")
