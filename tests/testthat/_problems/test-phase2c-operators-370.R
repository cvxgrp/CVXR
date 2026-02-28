# Extracted from test-phase2c-operators.R:370

# test -------------------------------------------------------------------------
lhs <- Constant(matrix(c(1, 2, 3, 4), 2, 2))
rhs <- Variable(c(2, 1))
m <- MulExpression(lhs, rhs)
