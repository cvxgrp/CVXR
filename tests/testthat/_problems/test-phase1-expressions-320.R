# Extracted from test-phase1-expressions.R:320

# test -------------------------------------------------------------------------
c1 <- Constant(5, name = "five")
expect_equal(expr_name(c1), "five")
