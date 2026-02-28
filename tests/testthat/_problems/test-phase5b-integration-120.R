# Extracted from test-phase5b-integration.R:120

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
p <- Problem(Maximize(sum_entries(x)), list(x <= 1))
chain <- Chain(reductions = list(FlipObjective()))
