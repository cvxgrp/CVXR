# Extracted from test-phase4-constraints.R:385

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(3L, 1L))
y <- Variable(c(3L, 1L))
z <- Variable(c(3L, 1L))
constr <- ExpCone(x, y, z)
expect_equal(constr_size(constr), 9L)
