## CVXPY 1.9 NLP quasi-Newton (L-BFGS) tests
## (nlp_tests/test_quasi_newton.py::TestQuasiNewton).
##
## Verify IPOPT's L-BFGS mode (`hessian_approximation = "limited-memory"`), where
## the Hessian is approximated rather than computed exactly. These tests solve and
## assert the optimum only (no DerivativeChecker -- L-BFGS does not use the exact
## Hessian). CVXR requires an initial value() for every variable; the option
## threads through psolve(..., solver = "IPOPT", hessian_approximation = ...).

## @cvxpy nlp_tests/test_quasi_newton.py::TestQuasiNewton::test_rosenbrock_lbfgs
test_that("L-BFGS: Rosenbrock", {
  x <- Variable(2, name = "x"); value(x) <- c(-1.2, 1)
  prob <- Problem(Minimize((1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2))
  .de_ipopt_solve(prob, hessian_approximation = "limited-memory")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_quasi_newton.py::TestQuasiNewton::test_hs071_lbfgs
test_that("L-BFGS: HS071", {
  x <- Variable(4, bounds = list(0, 6)); value(x) <- c(1, 5, 5, 1)
  obj <- Minimize(x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3])
  prob <- Problem(obj, list(x[1] * x[2] * x[3] * x[4] >= 25, sum_entries(x^2) == 40))
  .de_ipopt_solve(prob, hessian_approximation = "limited-memory")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)),
               c(0.75450865, 4.63936861, 3.78856881, 1.88513184), tolerance = 1e-3)
})

## @cvxpy nlp_tests/test_quasi_newton.py::TestQuasiNewton::test_portfolio_lbfgs
test_that("L-BFGS: portfolio (quad_form)", {
  r <- c(0.026002150277777, 0.008101316405671, 0.073715909491990)
  Q <- matrix(c(0.018641039983891, 0.003598532927677, 0.001309759253660,
                0.003598532927677, 0.006436938322676, 0.004887265158407,
                0.001309759253660, 0.004887265158407, 0.068682765454814), 3, 3, byrow = TRUE)
  x <- Variable(3); value(x) <- c(10, 10, 10)
  prob <- Problem(Minimize(quad_form(x, Q)),
                  list(sum_entries(x) <= 1000, t(r) %*% x >= 50, x >= 0))
  .de_ipopt_solve(prob, hessian_approximation = "limited-memory")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(497.045504, 0.0, 502.954496), tolerance = 1e-3)
})

## @cvxpy nlp_tests/test_quasi_newton.py::TestQuasiNewton::test_analytic_center_lbfgs
## R note: A is random under R's RNG; CVXPY asserts only OPTIMAL status.
test_that("L-BFGS: analytic center (log barrier)", {
  set.seed(0); m <- 50; n <- 4
  b <- rep(1, m)
  rand <- matrix(rnorm((m - 2 * n) * n), m - 2 * n, n)
  A <- rbind(rand, diag(n), -diag(n))
  x <- Variable(n); value(x) <- rep(0, n)
  prob <- Problem(Minimize(-sum_entries(log(b - A %*% x))))
  .de_ipopt_solve(prob, hessian_approximation = "limited-memory")
  expect_equal(status(prob), "optimal")
})

## @cvxpy nlp_tests/test_quasi_newton.py::TestQuasiNewton::test_exact_vs_lbfgs_solution_quality
test_that("L-BFGS: exact vs limited-memory agree on Rosenbrock", {
  x <- Variable(2, name = "x")
  obj <- Minimize((1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2)
  value(x) <- c(-1.2, 1)
  .de_ipopt_solve(Problem(obj), hessian_approximation = "exact")
  x_exact <- as.numeric(value(x))
  value(x) <- c(-1.2, 1)
  .de_ipopt_solve(Problem(obj), hessian_approximation = "limited-memory")
  x_lbfgs <- as.numeric(value(x))
  expect_equal(x_exact, x_lbfgs, tolerance = 1e-4)
  expect_equal(x_exact, c(1, 1), tolerance = 1e-4)
  expect_equal(x_lbfgs, c(1, 1), tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_quasi_newton.py::TestQuasiNewton::test_large_scale_lbfgs
## R note: Q, c are random under R's RNG; CVXPY asserts only OPTIMAL status.
test_that("L-BFGS: larger quadratic", {
  set.seed(42); n <- 100
  Q <- matrix(rnorm(n * n), n, n); Q <- t(Q) %*% Q / n
  cc <- rnorm(n) / n
  x <- Variable(n); value(x) <- rep(1 / n, n)
  prob <- Problem(Minimize(0.5 * quad_form(x, Q) + t(cc) %*% x),
                  list(sum_entries(x) == 1, x >= 0))
  .de_ipopt_solve(prob, hessian_approximation = "limited-memory")
  expect_equal(status(prob), "optimal")
})

## @cvxpy nlp_tests/test_quasi_newton.py::TestQuasiNewton::test_socp_lbfgs
test_that("L-BFGS: SOCP", {
  x <- Variable(3); y <- Variable(); value(x) <- c(0, 0, 0); value(y) <- 1
  prob <- Problem(Minimize(3 * x[1] + 2 * x[2] + x[3]),
                  list(p_norm(x, 2) <= y, x[1] + x[2] + 3 * x[3] >= 1.0, y <= 5))
  .de_ipopt_solve(prob, hessian_approximation = "limited-memory")
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), -13.548638814247532, tolerance = 1e-3)
})

## @cvxpy nlp_tests/test_quasi_newton.py::TestQuasiNewton::test_constrained_log_lbfgs
test_that("L-BFGS: constrained log -> uniform", {
  n <- 20
  x <- Variable(n, pos = TRUE); value(x) <- rep(1 / n, n)
  prob <- Problem(Minimize(-sum_entries(log(x))), list(sum_entries(x) == 1))
  .de_ipopt_solve(prob, hessian_approximation = "limited-memory")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), rep(1 / n, n), tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_quasi_newton.py::TestQuasiNewton::test_entropy_lbfgs
## Nonconvex entropy minimization -> distribution concentrates on one point.
test_that("L-BFGS: entropy minimization (nonconvex)", {
  set.seed(0); n <- 10
  q <- Variable(n, nonneg = TRUE)
  qv <- runif(n); value(q) <- qv / sum(qv)
  prob <- Problem(Minimize(sum_entries(entr(q))), list(sum_entries(q) == 1))
  .de_ipopt_solve(prob, hessian_approximation = "limited-memory")
  expect_equal(sum(as.numeric(value(q)) > 1e-8), 1L)
})
