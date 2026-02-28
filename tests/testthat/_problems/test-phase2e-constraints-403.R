# Extracted from test-phase2e-constraints.R:403

# test -------------------------------------------------------------------------
x <- Variable(3)
y <- Variable(3)
constr <- Equality(x, y)
expect_error(constr_expr(constr), "ambiguous")
