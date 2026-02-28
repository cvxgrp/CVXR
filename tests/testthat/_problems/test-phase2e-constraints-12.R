# Extracted from test-phase2e-constraints.R:12

# test -------------------------------------------------------------------------
x <- Variable(3, name = "x")
constr <- Zero(x)
expect_true(S7_inherits(constr, Zero))
