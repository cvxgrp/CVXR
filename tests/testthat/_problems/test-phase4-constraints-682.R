# Extracted from test-phase4-constraints.R:682

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
W <- Variable(c(3L, 2L))
z <- Variable(c(2L, 1L))
alpha <- Constant(matrix(c(0.2, 0.3, 0.5, 0.4, 0.4, 0.2), 3, 2))
constr <- PowConeND(W, z, alpha, axis = 2L)
expect_equal(cone_sizes(constr), c(4L, 4L))
