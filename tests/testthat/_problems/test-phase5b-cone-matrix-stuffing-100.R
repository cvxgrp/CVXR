# Extracted from test-phase5b-cone-matrix-stuffing.R:100

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
p <- Problem(Minimize(2 * x), list(x >= 0, x <= 5))
p_dcp <- reduction_apply(Dcp2Cone(), p)[[1]]
