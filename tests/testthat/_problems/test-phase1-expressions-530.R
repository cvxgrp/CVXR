# Extracted from test-phase1-expressions.R:530

# test -------------------------------------------------------------------------
x <- Variable(c(3, 4))
expect_equal(expr_ndim(x), 2L)
