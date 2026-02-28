# Extracted from test-phase4-constraints.R:195

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
t_var <- Variable(1L)
x_var <- Variable(c(2L, 1L))
constr <- SOC(t_var, x_var)
save_dual_value(constr, c(10, 20, 30))
