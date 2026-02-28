# Extracted from test-phase4-constraints.R:83

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
t_var <- Variable(c(3L, 1L))
X_var <- Variable(c(3L, 2L))
constr <- SOC(t_var, X_var, axis = 1L)
expect_equal(constr@axis, 1L)
expect_equal(num_cones(constr), 3L)
