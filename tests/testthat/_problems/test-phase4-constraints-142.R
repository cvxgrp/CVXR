# Extracted from test-phase4-constraints.R:142

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
t_var <- Variable(1L)
x_var <- Variable(c(3L, 1L))
constr <- SOC(t_var, x_var)
expect_equal(cone_sizes(constr), 4L)
