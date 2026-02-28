# Extracted from test-phase4-constraints.R:738

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
t_var <- Variable(1L)
x_var <- Variable(c(3L, 1L))
soc <- SOC(t_var, x_var)
psd <- PSD(Variable(c(2L, 2L)))
ec <- ExpCone(Variable(1L), Variable(1L), Variable(1L))
pc3 <- PowCone3D(Variable(1L), Variable(1L), Variable(1L), 0.5)
for (c in list(soc, psd, ec, pc3)) {
    expect_true(S7::S7_inherits(c, Cone))
    expect_true(S7::S7_inherits(c, Constraint))
  }
