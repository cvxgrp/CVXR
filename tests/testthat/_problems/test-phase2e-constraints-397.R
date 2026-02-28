# Extracted from test-phase2e-constraints.R:397

# test -------------------------------------------------------------------------
x <- Variable(3, name = "x")
constr <- Zero(x)
expect_identical(constr_expr(constr)@id, x@id)
