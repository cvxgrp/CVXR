# Extracted from test-phase4-constraints.R:117

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
t_var <- Variable(1L)
x_var <- Variable(c(3L, 1L))
constr <- SOC(t_var, x_var)
nm <- expr_name(constr)
