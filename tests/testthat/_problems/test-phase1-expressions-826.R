# Extracted from test-phase1-expressions.R:826

# test -------------------------------------------------------------------------
x <- Variable(c(2, 2), symmetric = TRUE)
expect_true(is_hermitian(x))
