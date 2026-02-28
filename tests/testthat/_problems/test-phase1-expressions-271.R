# Extracted from test-phase1-expressions.R:271

# test -------------------------------------------------------------------------
c1 <- Constant(5)
expect_length(grad(c1), 0)
