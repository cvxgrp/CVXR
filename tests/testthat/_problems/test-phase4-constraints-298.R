# Extracted from test-phase4-constraints.R:298

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
X <- Variable(c(3L, 3L))
constr <- PSD(X)
dc <- dual_cone(constr)
