# Extracted from test-phase4-constraints.R:125

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
t_var <- Variable(1L)
x_var <- Variable(c(3L, 1L))
constr <- SOC(t_var, x_var, axis = 2L)
d <- get_data(constr)
