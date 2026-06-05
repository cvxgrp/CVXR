## CVXPY 1.9 NLP diff-engine stress tests for nonlinear compositions
## (nlp_tests/stress_tests_diff_engine/test_compositions.py).
##
## Each CVXPY test runs DerivativeChecker.run_and_assert (and some also solve +
## assert an optimum). The R port runs the faithful mirror `nlp_derivative_check`
## at an in-domain point set via value(); where CVXPY additionally asserts a
## solved optimum we add a guarded solve. cp.multiply -> `*`, cp.square -> `^2`,
## cp.matmul -> `%*%`, cp.nlp.sin/cos -> sin/cos, transpose -> t().

## @cvxpy nlp_tests/stress_tests_diff_engine/test_compositions.py::TestCompositions::test_left_matmul_composition
test_that("DE composition: Trace(exp(A @ X)) with linear eq", {
  set.seed(0); X <- Variable(c(10L, 10L), bounds = list(-0.2, 0.2))
  A <- matrix(runif(100), 10, 10)
  value(X) <- matrix(runif(100) * 0.4 - 0.2, 10, 10)
  ## CVXPY X[1,1]+X[2,2] -> R X[2,2]+X[3,3]
  prob <- Problem(Minimize(matrix_trace(exp(A %*% X))), list(X[2, 2] + X[3, 3] == 0.1))
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_compositions.py::TestCompositions::test_right_matmul_composition
test_that("DE composition: sum(exp(X @ A)) with linear eq", {
  set.seed(0); X <- Variable(c(10L, 10L), bounds = list(-0.2, 0.2))
  A <- matrix(runif(100), 10, 10)
  value(X) <- matrix(runif(100) * 0.4 - 0.2, 10, 10)
  prob <- Problem(Minimize(sum_entries(exp(X %*% A))), list(X[2, 2] + X[3, 3] == 0.1))
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_compositions.py::TestCompositions::test_multiply_linear_composition
test_that("DE composition: sum((A @ X) * (B @ Y))", {
  set.seed(0); m <- 20; n <- 5
  A <- matrix(runif(m * n), m, n); B <- matrix(runif(m * n), m, n)
  X <- Variable(c(n, n), bounds = list(-1, 1)); Y <- Variable(c(n, n), bounds = list(-1, 1))
  value(X) <- matrix(runif(n * n), n, n); value(Y) <- matrix(runif(n * n), n, n)
  prob <- Problem(Minimize(sum_entries((A %*% X) * (B %*% Y))))
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_compositions.py::TestCompositions::test_multiply_nonlinear_composition
test_that("DE composition: sum((A @ X)^2 * logistic(B @ Y))", {
  set.seed(0); m <- 20; n <- 5
  A <- matrix(runif(m * n), m, n); B <- matrix(runif(m * n), m, n)
  X <- Variable(c(n, n), bounds = list(-1, 1)); Y <- Variable(c(n, n), bounds = list(-1, 1))
  value(X) <- matrix(runif(n * n), n, n); value(Y) <- matrix(runif(n * n), n, n)
  prob <- Problem(Minimize(sum_entries((A %*% X)^2 * logistic(B %*% Y))))
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_compositions.py::TestCompositions::test_multiply_nonlinear_composition_transpose
test_that("DE composition: sum(t(A @ X)^2 * logistic(B @ Y))", {
  set.seed(0); m <- 10; n <- 10
  A <- matrix(runif(m * n), m, n); B <- matrix(runif(m * n), m, n)
  X <- Variable(c(n, n), bounds = list(-1, 1)); Y <- Variable(c(n, n), bounds = list(-1, 1))
  value(X) <- matrix(runif(n * n), n, n); value(Y) <- matrix(runif(n * n), n, n)
  prob <- Problem(Minimize(sum_entries(t(A %*% X)^2 * logistic(B %*% Y))))
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_compositions.py::TestCompositions::test_quad_form_composition
test_that("DE composition: quad_form(sin(x) * x, Q) indefinite Q", {
  set.seed(0); n <- 25
  Q <- matrix(runif(n * n), n, n); Q <- Q + t(Q)
  x <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n)
  nlp_derivative_check(Problem(Minimize(quad_form(sin(x) * x, Q))))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_compositions.py::TestCompositions::test_quad_form_composition_two
test_that("DE composition: quad_form(sin(x) * (x * y), Q)", {
  set.seed(0); n <- 10
  Q <- matrix(runif(n * n), n, n); Q <- Q + t(Q)
  x <- Variable(n, bounds = list(-1, 1)); y <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n); value(y) <- runif(n)
  nlp_derivative_check(Problem(Minimize(quad_form(sin(x) * (x * y), Q))))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_compositions.py::TestCompositions::test_matmul_composition_one
test_that("DE composition: sum(X %*% cos(Y))", {
  set.seed(0); m <- 5; n <- 7; p <- 11
  X <- Variable(c(m, n), bounds = list(-1, 1)); Y <- Variable(c(n, p), bounds = list(-2, 2))
  value(X) <- matrix(runif(m * n), m, n); value(Y) <- matrix(runif(n * p), n, p)
  nlp_derivative_check(Problem(Minimize(sum_entries(X %*% cos(Y)))))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_compositions.py::TestCompositions::test_matmul_composition_two
test_that("DE composition: sum((X %*% X) %*% (cos(Y) + X))", {
  set.seed(0); m <- 5; n <- 5; p <- 5
  X <- Variable(c(m, n), bounds = list(-1, 1)); Y <- Variable(c(n, p), bounds = list(-2, 2))
  value(X) <- matrix(runif(m * n), m, n); value(Y) <- matrix(runif(n * p), n, p)
  nlp_derivative_check(Problem(Minimize(sum_entries((X %*% X) %*% (cos(Y) + X)))))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_compositions.py::TestCompositions::test_matmul_composition_three
test_that("DE composition: sum((X %*% t(X)) %*% t(cos(Y) + X))", {
  set.seed(0); m <- 5; n <- 5; p <- 5
  X <- Variable(c(m, n), bounds = list(-1, 1)); Y <- Variable(c(n, p), bounds = list(-2, 2))
  value(X) <- matrix(runif(m * n), m, n); value(Y) <- matrix(runif(n * p), n, p)
  nlp_derivative_check(Problem(Minimize(sum_entries((X %*% t(X)) %*% t(cos(Y) + X)))))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_compositions.py::TestCompositions::test_div
## Python `1 / (x / A @ x)` parses as `1 / ((x / A) @ x)` (`/` and `@` share
## precedence, left-assoc); `x / A` broadcasts (n,1)/(n,n) -> (n,n). Maximizing
## the denominator drives x to its upper bound 2.
test_that("DE composition: sum(1 / ((x / A) @ x)) -> x = 2", {
  set.seed(0); n <- 5
  x <- Variable(c(n, 1L), bounds = list(1, 2))
  A <- matrix(runif(n * n), n, n)
  value(x) <- matrix(runif(n) + 1, n, 1)
  prob <- Problem(Minimize(sum_entries(1 / ((x / A) %*% x))))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  expect_equal(as.numeric(value(x)), rep(2, n), tolerance = 1e-3)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_compositions.py::TestCompositions::test_negative_power
test_that("DE composition: sum(power(A @ x, -0.5))", {
  set.seed(0); n <- 5
  x <- Variable(c(n, 1L), bounds = list(1, 2))
  A <- matrix(runif(n * n), n, n)
  value(x) <- matrix(runif(n) + 1, n, 1)
  prob <- Problem(Minimize(sum_entries(power(A %*% x, -0.5))))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  expect_equal(status(prob), "optimal")
})
