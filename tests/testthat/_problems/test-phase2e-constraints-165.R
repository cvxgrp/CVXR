# Extracted from test-phase2e-constraints.R:165

# test -------------------------------------------------------------------------
x <- Variable(3, name = "x")
constr <- NonPos(x)
expect_true(S7_inherits(constr, NonPos))
