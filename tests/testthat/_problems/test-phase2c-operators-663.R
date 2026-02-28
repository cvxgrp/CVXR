# Extracted from test-phase2c-operators.R:663

# test -------------------------------------------------------------------------
x <- Variable(c(3, 1))
c1 <- Constant(matrix(1:3, ncol = 1))
m <- Multiply(c1, x)
