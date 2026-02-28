# Extracted from test-phase2c-operators.R:387

# test -------------------------------------------------------------------------
lhs <- Constant(matrix(1:4, 2, 2))
rhs <- Constant(matrix(c(1, 0, 0, 1), 2, 2))
m <- MulExpression(lhs, rhs)
