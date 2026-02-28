# Extracted from test-phase4-constraints.R:542

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(2L, 1L))
y <- Variable(c(2L, 1L))
z <- Variable(c(2L, 1L))
constr <- PowCone3D(x, y, z, 0.5)
expect_equal(cone_sizes(constr), c(3L, 3L))
