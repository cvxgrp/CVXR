# Extracted from test-phase5b-cone-matrix-stuffing.R:32

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
z <- Zero(x - Constant(c(1, 2)))
nn <- NonNeg(x)
cmap <- group_constraints(list(z, nn))
