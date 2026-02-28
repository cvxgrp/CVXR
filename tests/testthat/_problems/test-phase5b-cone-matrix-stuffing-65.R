# Extracted from test-phase5b-cone-matrix-stuffing.R:65

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
p <- Problem(Minimize(sum_entries(x)), list(x >= 0))
cms <- ConeMatrixStuffing()
