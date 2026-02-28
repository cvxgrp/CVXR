# Extracted from test-phase2d-dispatch.R:356

# test -------------------------------------------------------------------------
x <- Variable(3, name = "x")
y <- Variable(3, name = "y")
expect_equal(expr_name(x + y), "x + y")
