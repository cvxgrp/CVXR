# Extracted from test-phase2e-constraints.R:351

# test -------------------------------------------------------------------------
x <- Variable(c(2, 1), name = "x")
z <- Variable(c(2, 1), name = "z")
constr <- Inequality(x, z)
A <- Variable(c(2, 2), name = "A")
B <- Variable(c(2, 2), name = "B")
copy <- expr_copy(constr, args = list(A, B))
