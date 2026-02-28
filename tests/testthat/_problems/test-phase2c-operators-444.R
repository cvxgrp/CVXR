# Extracted from test-phase2c-operators.R:444

# test -------------------------------------------------------------------------
v1 <- matrix(c(1, 2, 3), ncol = 1)
v2 <- matrix(c(10, 20, 30), ncol = 1)
m <- Multiply(Constant(v1), Constant(v2))
