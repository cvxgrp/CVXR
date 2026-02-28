# Extracted from test-phase1-expressions.R:596

# test -------------------------------------------------------------------------
x <- Variable(c(2, 2), symmetric = TRUE)
val <- matrix(c(1, 2, 3, 4), 2, 2)
result <- project(x, val)
