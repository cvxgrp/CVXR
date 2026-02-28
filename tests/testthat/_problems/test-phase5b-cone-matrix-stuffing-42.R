# Extracted from test-phase5b-cone-matrix-stuffing.R:42

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
t_var <- Variable(1)
x <- Variable(3)
soc <- SOC(t_var, x)
cmap <- group_constraints(list(soc))
