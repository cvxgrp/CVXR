# Extracted from test-phase5b-cone-matrix-stuffing.R:132

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
xv <- Variable(1)
yz <- Variable(2)
p <- Problem(Minimize(xv), list(SOC(xv, yz)))
p_dcp <- reduction_apply(Dcp2Cone(), p)[[1]]
