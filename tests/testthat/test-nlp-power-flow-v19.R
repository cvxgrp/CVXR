## CVXPY 1.9 NLP AC power-flow example (nlp_tests/test_power_flow.py).
##
## Dense AC optimal power flow (DNLP paper). All data is deterministic, so the
## minimum generation cost (3087.84) is reproducible verbatim. C = cos(theta -
## theta^T), S = sin(theta - theta^T) and the v v^T outer product make this a
## genuinely nonconvex NLP. Solve with IPOPT, assert the cost, and run the
## DerivativeChecker mirror. cp.multiply -> `*`, cp.sum(., axis=1) -> axis 1.

## @cvxpy nlp_tests/test_power_flow.py::TestPowerFlowIPOPT::test_power_flow_dense_formulation
test_that("power flow: dense AC OPF reaches the known optimum", {
  if (!.de_have_backend()) skip("No R NLP backend")
  skip_if_not_installed("ipopt")
  N <- 9
  p_min <- rep(0, N); p_max <- rep(0, N); q_min <- rep(0, N); q_max <- rep(0, N)
  p_min[c(1, 2, 3)] <- c(10, 10, 10); p_max[c(1, 2, 3)] <- c(250, 300, 270)
  q_min[c(1, 2, 3)] <- c(-5, -5, -5)
  p_min[c(5, 7, 9)] <- c(-54, -60, -75); p_max[c(5, 7, 9)] <- c(-54, -60, -75)
  q_min[c(5, 7, 9)] <- c(-18, -21, -30); q_max[c(5, 7, 9)] <- c(-18, -21, -30)
  v_min <- 0.9; v_max <- 1.1

  ## Branch data: (from_bus, to_bus, resistance, reactance, susceptance), 0-based buses.
  branch <- matrix(c(
    0, 3, 0.0,    0.0576, 0.0,
    3, 4, 0.017,  0.092,  0.158,
    5, 4, 0.039,  0.17,   0.358,
    2, 5, 0.0,    0.0586, 0.0,
    5, 6, 0.0119, 0.1008, 0.209,
    7, 6, 0.0085, 0.072,  0.149,
    1, 7, 0.0,    0.0625, 0.0,
    7, 8, 0.032,  0.161,  0.306,
    3, 8, 0.01,   0.085,  0.176), ncol = 5, byrow = TRUE)
  M <- nrow(branch); base_MVA <- 100
  from_bus <- branch[, 1] + 1L; to_bus <- branch[, 2] + 1L

  A <- matrix(0, N, M)
  for (k in seq_len(M)) { A[from_bus[k], k] <- 1; A[to_bus[k], k] <- -1 }
  z <- (branch[, 3] + 1i * branch[, 4]) / base_MVA
  Y_0 <- A %*% diag(1 / z) %*% t(A)
  y_sh <- 0.5 * (1i * branch[, 5]) * base_MVA
  Y_sh <- diag(diag(A %*% diag(y_sh) %*% t(A)))
  G <- Re(Y_0) + Re(Y_sh); B <- Im(Y_0) + Im(Y_sh)

  theta <- Variable(c(N, 1L)); P <- Variable(c(N, N)); Q <- Variable(c(N, N))
  v <- Variable(c(N, 1L), bounds = list(v_min, v_max))
  p <- Variable(N, bounds = list(p_min, p_max))
  q <- Variable(N, bounds = list(q_min, q_max))
  C <- cos(theta - t(theta)); S <- sin(theta - t(theta))
  vvt <- v %*% t(v)
  constr <- list(
    theta[1] == 0,
    p == sum_entries(P, axis = 1),
    q == sum_entries(Q, axis = 1),
    P == vvt * (G * C + B * S),
    Q == vvt * (G * S - B * C))
  cost <- 0.11 * p[1]^2 + 5 * p[1] + 150 + 0.085 * p[2]^2 + 1.2 * p[2] + 600 +
    0.1225 * p[3]^2 + p[3] + 335
  prob <- Problem(Minimize(cost), constr)
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_lt(abs(value(prob) - 3087.84) / value(prob), 1e-4)
  nlp_derivative_check(prob)
})
