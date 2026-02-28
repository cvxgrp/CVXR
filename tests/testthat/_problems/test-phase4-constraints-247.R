# Extracted from test-phase4-constraints.R:247

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
X <- Variable(c(3L, 3L))
constr <- PSD(X)
expect_equal(cone_sizes(constr), 3L)
