# Extracted from test-phase4-constraints.R:434

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(1L)
y <- Variable(1L)
z <- Variable(1L)
constr <- ExpCone(x, y, z)
dc <- dual_cone(constr)
