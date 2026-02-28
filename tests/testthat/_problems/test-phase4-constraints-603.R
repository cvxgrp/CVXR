# Extracted from test-phase4-constraints.R:603

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(1L)
y <- Variable(1L)
z <- Variable(1L)
constr <- PowCone3D(x, y, z, 0.5)
nm <- expr_name(constr)
