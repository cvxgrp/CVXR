# Extracted from test-phase4-constraints.R:75

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
t_var <- Variable(c(2L, 1L))
X_var <- Variable(c(3L, 2L))
constr <- SOC(t_var, X_var, axis = 2L)
expect_true(S7::S7_inherits(constr, SOC))
expect_equal(num_cones(constr), 2L)
