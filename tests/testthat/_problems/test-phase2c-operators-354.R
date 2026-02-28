# Extracted from test-phase2c-operators.R:354

# test -------------------------------------------------------------------------
lhs <- Constant(matrix(1:6, 3, 2))
rhs <- Variable(c(2, 1))
m <- MulExpression(lhs, rhs)
