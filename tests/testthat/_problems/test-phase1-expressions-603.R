# Extracted from test-phase1-expressions.R:603

# test -------------------------------------------------------------------------
x <- Variable(c(2, 2), PSD = TRUE)
val <- matrix(c(1, 0, 0, -1), 2, 2)
result <- project(x, val)
