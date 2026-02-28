# Extracted from test-phase2d-dispatch.R:361

# test -------------------------------------------------------------------------
x <- Variable(3, name = "x")
expect_equal(expr_name(-x), "-x")
