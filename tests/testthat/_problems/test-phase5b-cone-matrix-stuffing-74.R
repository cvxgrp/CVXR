# Extracted from test-phase5b-cone-matrix-stuffing.R:74

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
p <- Problem(Minimize(exp(x)))
cms <- ConeMatrixStuffing()
