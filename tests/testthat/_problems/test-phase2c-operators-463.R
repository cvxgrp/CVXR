# Extracted from test-phase2c-operators.R:463

# test -------------------------------------------------------------------------
lhs <- Variable(1, nonneg = TRUE)
rhs <- Variable(1, nonpos = TRUE)
m <- Multiply(lhs, rhs)
