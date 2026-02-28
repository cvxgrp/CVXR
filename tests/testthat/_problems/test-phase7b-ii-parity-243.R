# Extracted from test-phase7b-ii-parity.R:243

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
P <- -matrix(c(2, 0.5, 0.5, 1), 2, 2)
result <- decomp_quad(P)
