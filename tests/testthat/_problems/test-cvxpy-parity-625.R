# Extracted from test-cvxpy-parity.R:625

# prequel ----------------------------------------------------------------------
library(testthat)
library(CVXR)

# test -------------------------------------------------------------------------
n <- 3L
P_mat <- diag(n) + 1e-12 * matrix(rnorm(n * n), n, n)
P_mat <- (P_mat + t(P_mat)) / 2
P <- Constant(P_mat)
x <- Variable(n)
qf <- QuadForm(x, P)
