# Extracted from test-phase2c-operators.R:335

# test -------------------------------------------------------------------------
lhs <- Variable(c(3, 4))
rhs <- Constant(matrix(1:8, 4, 2))
m <- MulExpression(lhs, rhs)
