# Extracted from test-phase4-constraints.R:369

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(3L, 1L))
y <- Variable(c(3L, 1L))
z <- Variable(c(3L, 1L))
constr <- ExpCone(x, y, z)
expect_equal(num_cones(constr), 3L)
