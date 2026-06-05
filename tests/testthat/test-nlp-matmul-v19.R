## CVXPY 1.9 DNLP matmul parity (nlp_tests/test_matmul.py).
##
## Matrix multiplication in a DNLP: variable-times-variable (X %*% Y),
## constant-factor (const %*% f(Y), f(X) %*% const), and 1-D operands. All are
## status-only solves in CVXPY plus a Python DerivativeChecker (dropped here).
## Variables use box bounds (CVXR `bounds = list(lo, hi)`); cp.nlp.sin/cos map
## to CVXR sin/cos.

.mm_nlp_optimal <- function(prob) {
  skip_if_not_installed("sparsediff")
  if (!requireNamespace("ipopt", quietly = TRUE) &&
      !requireNamespace("Uno", quietly = TRUE)) {
    skip("No R NLP backend (ipopt or Uno) installed")
  }
  output <- capture.output(psolve(prob, nlp = TRUE), type = "output")
  expect_equal(status(prob), "optimal")
}

## @cvxpy nlp_tests/test_matmul.py::TestMatmul::test_simple_matmul_graph_form
test_that("DNLP matmul: t == sum(X %*% Y) graph form solves", {
  set.seed(0)
  m <- 5; n <- 7; p <- 11
  X <- Variable(c(m, n), bounds = list(-1, 1))
  Y <- Variable(c(n, p), bounds = list(-2, 2))
  tt <- Variable()
  value(X) <- matrix(runif(m * n), m, n)
  value(Y) <- matrix(runif(n * p), n, p)
  value(tt) <- sum(value(X) %*% value(Y))
  .mm_nlp_optimal(Problem(Minimize(tt), list(tt == sum_entries(X %*% Y))))
})

## @cvxpy nlp_tests/test_matmul.py::TestMatmul::test_simple_matmul_not_graph_form
test_that("DNLP matmul: sum(X %*% Y) not-graph form solves", {
  set.seed(0)
  m <- 5; n <- 7; p <- 11
  X <- Variable(c(m, n), bounds = list(-1, 1))
  Y <- Variable(c(n, p), bounds = list(-2, 2))
  value(X) <- matrix(runif(m * n), m, n)
  value(Y) <- matrix(runif(n * p), n, p)
  .mm_nlp_optimal(Problem(Minimize(sum_entries(X %*% Y))))
})

## @cvxpy nlp_tests/test_matmul.py::TestMatmul::test_matmul_with_function_right
test_that("DNLP matmul: const %*% cos(Y) solves", {
  set.seed(0)
  m <- 5; n <- 7; p <- 11
  A <- matrix(runif(m * n), m, n)
  Y <- Variable(c(n, p), bounds = list(-2, 2))
  value(Y) <- matrix(runif(n * p), n, p)
  .mm_nlp_optimal(Problem(Minimize(sum_entries(A %*% cos(Y)))))
})

## @cvxpy nlp_tests/test_matmul.py::TestMatmul::test_matmul_with_function_left
test_that("DNLP matmul: cos(X) %*% const solves", {
  set.seed(0)
  m <- 5; n <- 7; p <- 11
  X <- Variable(c(m, n), bounds = list(-2, 2))
  B <- matrix(runif(n * p), n, p)
  value(X) <- matrix(runif(m * n), m, n)
  .mm_nlp_optimal(Problem(Minimize(sum_entries(cos(X) %*% B))))
})

## @cvxpy nlp_tests/test_matmul.py::TestMatmul::test_matmul_with_functions_both_sides
test_that("DNLP matmul: cos(X) %*% sin(Y) solves", {
  set.seed(0)
  m <- 5; n <- 7; p <- 11
  X <- Variable(c(m, n), bounds = list(-2, 2))
  Y <- Variable(c(n, p), bounds = list(-2, 2))
  value(X) <- matrix(runif(m * n), m, n)
  value(Y) <- matrix(runif(n * p), n, p)
  .mm_nlp_optimal(Problem(Minimize(sum_entries(cos(X) %*% sin(Y)))))
})

## @cvxpy nlp_tests/test_matmul.py::TestMatmul::test_matmul_1d_left_constant
test_that("DNLP matmul: 1-D constant on the left (a %*% sin(X)) solves", {
  set.seed(0)
  n <- 4; p <- 5
  a <- runif(n)
  X <- Variable(c(n, p), bounds = list(-1, 1))
  value(X) <- matrix(runif(n * p), n, p)
  .mm_nlp_optimal(Problem(Minimize(sum_entries(t(a) %*% sin(X)))))
})

## @cvxpy nlp_tests/test_matmul.py::TestMatmul::test_matmul_1d_right_constant
test_that("DNLP matmul: 1-D constant on the right (sin(X) %*% b) solves", {
  set.seed(0)
  m <- 5; n <- 4
  b <- runif(n)
  X <- Variable(c(m, n), bounds = list(-1, 1))
  value(X) <- matrix(runif(m * n), m, n)
  .mm_nlp_optimal(Problem(Minimize(sum_entries(sin(X) %*% b))))
})

## @cvxpy nlp_tests/test_matmul.py::TestMatmul::test_matmul_1d_dot
test_that("DNLP matmul: 1-D dot product (a %*% sin(x)) solves", {
  set.seed(0)
  n <- 6
  a <- runif(n)
  x <- Variable(n, bounds = list(-1, 1))
  value(x) <- runif(n)
  .mm_nlp_optimal(Problem(Minimize(t(a) %*% sin(x))))
})

## @cvxpy nlp_tests/test_matmul.py::TestMatmul::test_matmul_param_inside_transpose
## A parametrized matrix factor (A.T %*% sin(x) with A a Parameter) routes
## through the dense param-aware matmul (sd_left_matmul_dense fed the
## parameter's diff-engine node), so a re-solve after mutating A.value updates
## the matmul. DEVIATION: CVXPY asserts the parametric and hardcoded solutions
## are bit-identical (norm == 0); the parametric path uses a different
## (dense-param) matmul node than the hardcoded (sparse-constant) path, so we
## assert agreement to a tight tolerance. CVXPY clears x.value before each
## solve; the R NLP path requires an initial point, so we set the same start
## before each solve.
test_that("DNLP matmul: parametrized matrix factor (A.T %*% sin(x)) re-solves", {
  .mm_nlp_solve <- function(prob, x0) {
    skip_if_not_installed("sparsediff")
    if (!requireNamespace("ipopt", quietly = TRUE) &&
        !requireNamespace("Uno", quietly = TRUE)) {
      skip("No R NLP backend (ipopt or Uno) installed")
    }
    value(x0[[1L]]) <- x0[[2L]]
    capture.output(psolve(prob, nlp = TRUE, solver = "IPOPT"), type = "output")
    expect_equal(status(prob), "optimal")
    as.numeric(value(x0[[1L]]))
  }
  set.seed(0)
  m <- 4; p <- 5
  A1 <- matrix(runif(m * p), m, p)
  A2 <- matrix(runif(m * p), m, p)
  x <- Variable(m, bounds = list(-1, 1))
  start <- runif(m)

  hardcoded_sol1 <- .mm_nlp_solve(Problem(Minimize(sum_entries(t(A1) %*% sin(x)))), list(x, start))
  hardcoded_sol2 <- .mm_nlp_solve(Problem(Minimize(sum_entries(t(A2) %*% sin(x)))), list(x, start))

  A <- Parameter(c(m, p))
  prob <- Problem(Minimize(sum_entries(t(A) %*% sin(x))))
  value(A) <- A1
  param_sol1 <- .mm_nlp_solve(prob, list(x, start))
  value(A) <- A2
  param_sol2 <- .mm_nlp_solve(prob, list(x, start))

  expect_equal(param_sol1, hardcoded_sol1, tolerance = 1e-6)
  expect_equal(param_sol2, hardcoded_sol2, tolerance = 1e-6)
})
