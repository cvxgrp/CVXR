# Extracted from test-cvxpy-parity-phase7b.R:110

# test -------------------------------------------------------------------------
A <- Variable(c(2, 2))
B <- Variable(c(2, 2))
constr <- PSD(A - B)
expect_equal(constr@shape, c(2L, 2L))
save_leaf_value(A, matrix(c(3, 0, 0, 3), 2, 2))
