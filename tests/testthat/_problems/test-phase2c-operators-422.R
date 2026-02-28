# Extracted from test-phase2c-operators.R:422

# test -------------------------------------------------------------------------
lhs <- Constant(matrix(1:3, ncol = 1))
rhs <- Variable(c(3, 1))
m <- Multiply(lhs, rhs)
