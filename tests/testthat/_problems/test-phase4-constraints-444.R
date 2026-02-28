# Extracted from test-phase4-constraints.R:444

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(1L)
y <- Variable(1L)
z <- Variable(1L)
constr <- ExpCone(x, y, z)
save_dual_value(constr, c(10, 20, 30))
