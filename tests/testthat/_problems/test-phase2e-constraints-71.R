# Extracted from test-phase2e-constraints.R:71

# test -------------------------------------------------------------------------
x <- Variable(c(2, 1), name = "x")
z <- Variable(c(2, 1), name = "z")
constr <- Equality(x, z)
expect_true(S7_inherits(constr, Equality))
