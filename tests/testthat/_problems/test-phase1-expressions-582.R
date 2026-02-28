# Extracted from test-phase1-expressions.R:582

# test -------------------------------------------------------------------------
x <- Variable(3, nonneg = TRUE)
val <- matrix(c(-1, 2, -3))
result <- project(x, val)
