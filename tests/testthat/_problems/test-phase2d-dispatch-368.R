# Extracted from test-phase2d-dispatch.R:368

# test -------------------------------------------------------------------------
x <- Variable(3, name = "x")
y <- Variable(3, name = "y")
result <- x - y
expect_equal(expr_name(result), "x + -y")
