# Extracted from test-phase1-expressions.R:315

# test -------------------------------------------------------------------------
m <- matrix(1:6, 2, 3)
c1 <- Constant(m)
expect_true(grepl("2x3", expr_name(c1)))
