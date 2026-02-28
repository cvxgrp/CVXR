# Extracted from test-phase4-constraints.R:534

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(2L, 1L))
y <- Variable(c(2L, 1L))
z <- Variable(c(2L, 1L))
constr <- PowCone3D(x, y, z, 0.5)
expect_equal(num_cones(constr), 2L)
