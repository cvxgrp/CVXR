# Extracted from test-phase4-constraints.R:550

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(c(2L, 1L))
y <- Variable(c(2L, 1L))
z <- Variable(c(2L, 1L))
constr <- PowCone3D(x, y, z, 0.5)
expect_equal(constr_size(constr), 6L)
