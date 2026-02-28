# Extracted from test-cvxpy-parity-phase7b.R:66

# test -------------------------------------------------------------------------
x <- Variable(2)
z <- Variable(2)
constr <- x <= z
save_leaf_value(x, matrix(c(1, 1), 2, 1))
