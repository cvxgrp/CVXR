# Extracted from test-phase2e-constraints.R:139

# test -------------------------------------------------------------------------
x <- Variable(c(2, 1), name = "x")
z <- Variable(c(2, 1), name = "z")
constr <- Equality(x, z)
copy <- expr_copy(constr)
expect_true(S7_inherits(copy, Equality))
