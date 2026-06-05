## CVXPY 1.9 NLP example problems (nlp_tests/test_nlp_solvers.py::TestNLPExamples).
##
## Classic nonlinear programs from the IPOPT / JuMP documentation. CVXPY
## parametrizes over IPOPT/KNITRO/UNO/COPT; CVXR solves with the IPOPT backend
## (the robust interior-point reference -- CVXPY skips several of these for UNO).
## Each test solves, asserts the optimum, and runs the faithful DerivativeChecker
## mirror `nlp_derivative_check` (helper-nlp-diff-engine.R). CVXR requires an
## initial value() for every variable (unlike CVXPY/IPOPT's auto-init), so we set
## one. Value-asserts that depend on numpy's RNG-generated data are replaced by a
## status + KKT (derivative) check with a documented-deviation note; deterministic
## and noiseless problems keep their exact value-asserts. circle_packing tests
## live in test-nlp-solvers-circle-packing-v19.R.

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_hs071
test_that("NLP example: HS071", {
  x <- Variable(4, bounds = list(0, 6)); value(x) <- c(1, 5, 5, 1)
  obj <- Minimize(x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3])
  prob <- Problem(obj, list(x[1] * x[2] * x[3] * x[4] >= 25, sum_entries(x^2) == 40))
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)),
               c(0.75450865, 4.63936861, 3.78856881, 1.88513184), tolerance = 1e-5)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_mle
## R note: data is rnorm(1000) under R's RNG, not numpy's, so the numpy optimum
## (sigma=0.7708, mu=0.5941) is not reproducible. We assert the feasibility
## relation mu = sigma^2 plus the derivative (KKT) check instead.
test_that("NLP example: maximum likelihood with mu = sigma^2", {
  n <- 1000; set.seed(1234); data <- rnorm(n)
  mu <- Variable(1, name = "mu"); value(mu) <- 0
  sigma <- Variable(1, name = "sigma"); value(sigma) <- 1
  log_lik <- (n / 2) * log(1 / (2 * pi * sigma^2)) -
    sum_entries(square(data - mu)) / (2 * sigma^2)
  prob <- Problem(Maximize(log_lik), list(mu == sigma^2))
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(mu)), as.numeric(value(sigma))^2, tolerance = 1e-5)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_portfolio_opt
test_that("NLP example: portfolio (quad_form, r @ x)", {
  r <- c(0.026002150277777, 0.008101316405671, 0.073715909491990)
  Q <- matrix(c(0.018641039983891, 0.003598532927677, 0.001309759253660,
                0.003598532927677, 0.006436938322676, 0.004887265158407,
                0.001309759253660, 0.004887265158407, 0.068682765454814), 3, 3, byrow = TRUE)
  x <- Variable(3); value(x) <- c(10, 10, 10)
  prob <- Problem(Minimize(quad_form(x, Q)),
                  list(sum_entries(x) <= 1000, t(r) %*% x >= 50, x >= 0))
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(497.045504, 0.0, 502.954496), tolerance = 1e-3)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_portfolio_opt_sum_multiply
test_that("NLP example: portfolio (quad_form, sum(r * x))", {
  r <- c(0.026002150277777, 0.008101316405671, 0.073715909491990)
  Q <- matrix(c(0.018641039983891, 0.003598532927677, 0.001309759253660,
                0.003598532927677, 0.006436938322676, 0.004887265158407,
                0.001309759253660, 0.004887265158407, 0.068682765454814), 3, 3, byrow = TRUE)
  x <- Variable(3); value(x) <- c(10, 10, 10)
  prob <- Problem(Minimize(quad_form(x, Q)),
                  list(sum_entries(x) <= 1000, sum_entries(r * x) >= 50, x >= 0))
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(497.045504, 0.0, 502.954496), tolerance = 1e-3)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_rosenbrock
test_that("NLP example: Rosenbrock", {
  x <- Variable(2, name = "x"); value(x) <- c(-1.2, 1)
  prob <- Problem(Minimize((1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2))
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-5)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_socp
test_that("NLP example: SOCP via norm constraint", {
  x <- Variable(3); y <- Variable(); value(x) <- c(0, 0, 0); value(y) <- 1
  prob <- Problem(Minimize(3 * x[1] + 2 * x[2] + x[3]),
                  list(p_norm(x, 2) <= y, x[1] + x[2] + 3 * x[3] >= 1.0, y <= 5))
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), -13.548638814247532, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(-3.87462191, -2.12978826, 2.33480343), tolerance = 1e-4)
  expect_equal(as.numeric(value(y)), 5, tolerance = 1e-4)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_analytic_polytope_center_x_column_vector
## R note: A is random under R's RNG; CVXPY asserts only OPTIMAL status (no value).
test_that("NLP example: analytic polytope center (column x)", {
  set.seed(0); m <- 50; n <- 4
  b <- matrix(1, m, 1)
  rand <- matrix(rnorm((m - 2 * n) * n), m - 2 * n, n)
  A <- rbind(rand, diag(n), -diag(n))
  x <- Variable(c(n, 1L)); value(x) <- matrix(0, n, 1)
  prob <- Problem(Minimize(-sum_entries(log(b - A %*% x))))
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_portfolio_socp
## R note: mu/Sigma are random under R's RNG, so CVXPY's value (-1.934) is not
## reproducible; assert OPTIMAL status + KKT (derivative) check.
test_that("NLP example: portfolio SOCP (vector x)", {
  set.seed(858); n <- 100
  mu <- rnorm(n)
  Sigma <- matrix(rnorm(n * n), n, n); Sigma <- t(Sigma) %*% Sigma
  gamma <- 0.1; L <- t(chol(Sigma))  # lower
  x <- Variable(n, name = "x"); tt <- Variable(name = "t", bounds = list(0, Inf))
  value(x) <- rep(1 / n, n); value(tt) <- 1
  prob <- Problem(Minimize(-t(mu) %*% x + gamma * tt),
                  list(p_norm(t(L) %*% x, 2) <= tt, sum_entries(x) == 1, x >= 0))
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_portfolio_socp_x_column_vector
test_that("NLP example: portfolio SOCP (column x, sum(mu * x))", {
  set.seed(858); n <- 100
  mu <- matrix(rnorm(n), n, 1)
  Sigma <- matrix(rnorm(n * n), n, n); Sigma <- t(Sigma) %*% Sigma
  gamma <- 0.1; L <- t(chol(Sigma))
  x <- Variable(c(n, 1L), name = "x"); tt <- Variable(name = "t", bounds = list(0, Inf))
  value(x) <- matrix(1 / n, n, 1); value(tt) <- 1
  prob <- Problem(Minimize(-sum_entries(mu * x) + gamma * tt),
                  list(p_norm(t(L) %*% x, 2) <= tt, sum_entries(x) == 1, x >= 0))
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_localization
## Noiseless: rho_i = ||a_i - x_true||, so x_true is the exact minimizer for any
## a -- the value-assert reproduces under R's RNG.
test_that("NLP example: sensor localization (vector x)", {
  set.seed(42); m <- 10; x_true <- c(2.0, -1.5)
  a <- matrix(runif(m * 2, -5, 5), m, 2)
  rho <- sqrt(rowSums((a - matrix(x_true, m, 2, byrow = TRUE))^2))
  x <- Variable(2, name = "x"); t <- Variable(m, name = "t")
  value(x) <- c(0, 0); value(t) <- sqrt(rowSums(a^2))
  con <- list(t == sqrt(sum_entries(square(t(x) - a), axis = 1)))
  prob <- Problem(Minimize(sum_squares(t - rho)), con)
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), x_true, tolerance = 1e-4)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_localization2
test_that("NLP example: sensor localization (row x)", {
  set.seed(42); m <- 10; x_true <- c(2.0, -1.5)
  a <- matrix(runif(m * 2, -5, 5), m, 2)
  rho <- sqrt(rowSums((a - matrix(x_true, m, 2, byrow = TRUE))^2))
  x <- Variable(c(1L, 2L), name = "x"); t <- Variable(m, name = "t")
  value(x) <- matrix(0, 1, 2); value(t) <- sqrt(rowSums(a^2))
  con <- list(t == sqrt(sum_entries(square(x - a), axis = 1)))
  prob <- Problem(Minimize(sum_squares(t - rho)), con)
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), x_true, tolerance = 1e-4)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_geo_mean
test_that("NLP example: geometric mean", {
  x <- Variable(3, nonneg = TRUE); value(x) <- rep(1 / 3, 3)
  prob <- Problem(Maximize(geo_mean(x)), list(sum_entries(x) == 1))
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), rep(1 / 3, 3), tolerance = 1e-4)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_geo_mean2
test_that("NLP example: weighted geometric mean", {
  p <- c(.07, .12, .23, .19, .39)
  x <- Variable(5, nonneg = TRUE); value(x) <- rep(0.2, 5)
  prob <- Problem(Maximize(geo_mean(x, p)), list(sum_entries(x) <= 1))
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), p / sum(p), tolerance = 1e-4)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_div_composition
test_that("NLP example: maximize exp(1 / x)", {
  x <- Variable(nonneg = TRUE, bounds = list(1, 5)); value(x) <- 2
  prob <- Problem(Maximize(exp(1 / x)))
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-4)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_clnlbeam
## CVXPY skips its DerivativeChecker here (>10s on N=1000); we do too.
test_that("NLP example: clnlbeam", {
  N <- 1000; h <- 1 / N; alpha <- 350
  t <- Variable(N + 1, bounds = list(-1, 1)); value(t) <- rep(0, N + 1)
  x <- Variable(N + 1, bounds = list(-0.05, 0.05)); value(x) <- rep(0, N + 1)
  u <- Variable(N + 1); value(u) <- rep(0, N + 1)
  hi <- 2:(N + 1); lo <- 1:N
  control_terms <- 0.5 * h * (u[hi]^2 + u[lo]^2)
  trig_terms <- 0.5 * alpha * h * (cos(t[hi]) + cos(t[lo]))
  obj <- Minimize(sum_entries(control_terms + trig_terms))
  con <- list(
    x[hi] - x[lo] - 0.5 * h * (sin(t[hi]) + sin(t[lo])) == 0,
    t[hi] - t[lo] - 0.5 * h * (u[hi] + u[lo]) == 0)
  prob <- Problem(obj, con)
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 350.0, tolerance = 1e-1)
})
