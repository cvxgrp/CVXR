## CVXPY 1.9 NLP Sharpe-ratio example (nlp_tests/test_Sharpe_ratio.py).
##
## Maximizing the (nonconvex) Sharpe ratio square(mu'x)/quad_form(x, Sigma) over
## the simplex and the convex reformulation (min variance s.t. mu'x == 1, then
## renormalize) reach the same Sharpe ratio. That equality is a data-independent
## invariant, so R-generated Sigma/mu verify the same property CVXPY checks.

## @cvxpy nlp_tests/test_Sharpe_ratio.py::TestSharpeRatio::test_formulation_one
test_that("Sharpe ratio: nonconvex NLP matches convex reformulation", {
  if (!.de_have_backend()) skip("No R NLP backend")
  skip_if_not_installed("ipopt")
  set.seed(0); n <- 100
  Sigma <- matrix(runif(n * n), n, n); Sigma <- Sigma %*% t(Sigma)
  mu <- runif(n)
  x <- Variable(n, nonneg = TRUE)
  value(x) <- rep(1 / n, n)   # robust initialization (CVXPY does the same)
  prob <- Problem(Maximize(square(t(mu) %*% x) / quad_form(x, Sigma)),
                  list(sum_entries(x) == 1))
  .de_ipopt_solve(prob, hessian_approximation = "exact")
  x_noncvx <- as.numeric(value(x))

  prob_cvx <- Problem(Minimize(quad_form(x, Sigma)), list(t(mu) %*% x == 1))
  psolve(prob_cvx, solver = "CLARABEL")
  x_cvx <- as.numeric(value(x)); x_cvx <- x_cvx / sum(x_cvx)

  sharpe_nlp <- sum(mu * x_noncvx) / sqrt(as.numeric(t(x_noncvx) %*% Sigma %*% x_noncvx))
  sharpe_cvx <- sum(mu * x_cvx) / sqrt(as.numeric(t(x_cvx) %*% Sigma %*% x_cvx))
  expect_lt(abs(sharpe_nlp - sharpe_cvx), 1e-6)
  value(x) <- x_noncvx
  nlp_derivative_check(prob)
})
