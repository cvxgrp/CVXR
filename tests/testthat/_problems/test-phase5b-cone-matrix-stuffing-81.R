# Extracted from test-phase5b-cone-matrix-stuffing.R:81

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
p_dcp <- reduction_apply(Dcp2Cone(), Problem(Minimize(sum_entries(x)), list(x >= 1)))[[1]]
