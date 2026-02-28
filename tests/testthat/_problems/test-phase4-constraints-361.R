# Extracted from test-phase4-constraints.R:361

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
x <- Variable(1L)
y <- Variable(1L)
z <- Variable(1L)
constr <- ExpCone(x, y, z)
expect_false(is_dgp(constr))
