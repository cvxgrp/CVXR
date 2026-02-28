# Extracted from test-phase1-expressions.R:346

# test -------------------------------------------------------------------------
p <- Parameter(c(2, 2), name = "alpha")
expect_equal(expr_name(p), "alpha")
