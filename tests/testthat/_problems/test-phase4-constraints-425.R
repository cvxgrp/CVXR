# Extracted from test-phase4-constraints.R:425

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(1L)
y <- Variable(1L)
z <- Variable(1L)
constr <- ExpCone(x, y, z)
nm <- expr_name(constr)
