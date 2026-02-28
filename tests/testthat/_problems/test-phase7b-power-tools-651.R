# Extracted from test-phase7b-power-tools.R:651

# test -------------------------------------------------------------------------
t_var <- Variable(c(1L, 1L))
x1 <- Variable(c(1L, 1L))
x2 <- Variable(c(1L, 1L))
constrs <- powcone_constrs(t_var, list(x1, x2), 0.5)
