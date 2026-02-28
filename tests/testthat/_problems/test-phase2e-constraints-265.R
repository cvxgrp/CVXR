# Extracted from test-phase2e-constraints.R:265

# test -------------------------------------------------------------------------
x <- Variable(c(2, 1), name = "x")
z <- Variable(c(2, 1), name = "z")
constr <- Inequality(x, z)
expect_true(S7_inherits(constr, Inequality))
