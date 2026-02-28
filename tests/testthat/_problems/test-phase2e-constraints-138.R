# Extracted from test-phase2e-constraints.R:138

# test -------------------------------------------------------------------------
x <- Variable(c(2, 1), name = "x")
z <- Variable(c(2, 1), name = "z")
constr <- Equality(x, z)
copy <- expr_copy(constr)
