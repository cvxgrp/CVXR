# Extracted from test-cvxpy-parity-phase7b.R:121

# test -------------------------------------------------------------------------
A <- Variable(c(2, 2))
B <- Variable(c(2, 2))
constr <- PSD(A - B)
save_leaf_value(A, matrix(c(0, 0, 0, 0), 2, 2))
