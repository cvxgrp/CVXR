# Extracted from test-phase5b-cone-matrix-stuffing.R:147

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
p <- Problem(Minimize(x + 10))
p_dcp <- reduction_apply(Dcp2Cone(), p)[[1]]
