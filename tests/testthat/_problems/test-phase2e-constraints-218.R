# Extracted from test-phase2e-constraints.R:218

# test -------------------------------------------------------------------------
x <- Variable(3, name = "x")
constr <- NonNeg(x)
expect_true(S7_inherits(constr, NonNeg))
