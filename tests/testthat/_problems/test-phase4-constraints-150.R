# Extracted from test-phase4-constraints.R:150

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
t_var <- Variable(c(2L, 1L))
X_var <- Variable(c(3L, 2L))
constr <- SOC(t_var, X_var, axis = 2L)
expect_equal(constr_size(constr), 8L)
