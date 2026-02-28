# Extracted from test-phase2e-constraints.R:365

# test -------------------------------------------------------------------------
x <- Variable(c(2, 1), name = "x")
z <- Variable(c(2, 1), name = "z")
constr <- Inequality(x, z)
expect_equal(expr_name(constr), "x <= z")
