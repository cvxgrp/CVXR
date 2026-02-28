# Extracted from test-phase2e-constraints.R:339

# test -------------------------------------------------------------------------
x <- Variable(c(2, 1), name = "x")
z <- Variable(c(2, 1), name = "z")
constr <- Inequality(x, z)
copy <- expr_copy(constr)
expect_true(S7_inherits(copy, Inequality))
