# Extracted from test-phase5b-cone-matrix-stuffing.R:116

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
xy <- Variable(2)
p <- Problem(Minimize(sum_entries(xy)), list(xy[1] + xy[2] == 5, xy >= 0))
p_dcp <- reduction_apply(Dcp2Cone(), p)[[1]]
