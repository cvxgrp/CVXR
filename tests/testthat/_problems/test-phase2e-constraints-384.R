# Extracted from test-phase2e-constraints.R:384

# test -------------------------------------------------------------------------
x <- Variable(3)
constr <- Zero(x)
expect_true(is_real(constr))
