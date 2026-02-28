# Extracted from test-phase4-constraints.R:312

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(2L, 1L))
y <- Variable(c(2L, 1L))
z <- Variable(c(2L, 1L))
constr <- ExpCone(x, y, z)
expect_true(S7::S7_inherits(constr, ExpCone))
expect_true(S7::S7_inherits(constr, Cone))
