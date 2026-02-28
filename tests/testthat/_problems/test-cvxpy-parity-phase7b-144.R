# Extracted from test-cvxpy-parity-phase7b.R:144

# test -------------------------------------------------------------------------
x <- Variable(2)
z <- Variable(2)
a <- Variable(1)
b <- Variable(1)
exp_val <- x + z
scalar_exp <- a + b
constr <- SOC(scalar_exp, exp_val)
expect_equal(cone_sizes(constr), 3L)
