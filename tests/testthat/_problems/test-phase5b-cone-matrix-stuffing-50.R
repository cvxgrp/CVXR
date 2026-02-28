# Extracted from test-phase5b-cone-matrix-stuffing.R:50

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
nn <- NonNeg(x)
cmap <- group_constraints(list(nn))
