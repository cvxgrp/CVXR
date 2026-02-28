# Extracted from test-phase4-constraints.R:463

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(2L, 1L))
y <- Variable(c(2L, 1L))
z <- Variable(c(2L, 1L))
constr <- PowCone3D(x, y, z, 0.5)
expect_true(S7::S7_inherits(constr, PowCone3D))
expect_true(S7::S7_inherits(constr, Cone))
