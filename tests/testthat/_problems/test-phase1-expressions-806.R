# Extracted from test-phase1-expressions.R:806

# test -------------------------------------------------------------------------
x <- Variable(3, bounds = c(-1, 2))
val <- matrix(c(-5, 0, 10))
result <- project(x, val)
