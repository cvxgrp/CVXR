## CVXPY 1.9 NLP diff-engine stress tests for affine *vector* atoms
## (nlp_tests/stress_tests_diff_engine/test_affine_vector_atoms.py).
##
## Each CVXPY test builds a problem, runs a Python DerivativeChecker (analytic
## vs finite-difference objective gradient / constraint Jacobian / Lagrangian
## Hessian), then solves it as an NLP and asserts the optimum. The R port runs
## the faithful DerivativeChecker mirror `nlp_derivative_check()` at an in-domain
## point (helper-nlp-diff-engine.R) and, where R's shape semantics match CVXPY's,
## also solves and asserts the optimum. Box bounds use `bounds = list(lo, hi)`.
## `.de_nlp_solve` / `nlp_derivative_check` live in helper-nlp-diff-engine.R.

## @cvxpy nlp_tests/stress_tests_diff_engine/test_affine_vector_atoms.py::TestAffineDiffEngine::test_row_broadcast
test_that("DE affine-vector: row broadcast sum(x + Y)", {
  set.seed(0); m <- 3; n <- 4
  x <- Variable(c(1L, n), bounds = list(-2, 2))
  Y <- Variable(c(m, n), bounds = list(-1, 1))
  value(x) <- matrix(runif(n), 1, n)
  value(Y) <- matrix(runif(m * n), m, n)
  prob <- Problem(Minimize(sum_entries(x + Y)))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_true(all(abs(as.numeric(value(x)) - (-2)) < 1e-4))
  expect_true(all(abs(as.numeric(value(Y)) - (-1)) < 1e-4))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_affine_vector_atoms.py::TestAffineDiffEngine::test_col_broadcast
test_that("DE affine-vector: column broadcast sum(x + Y)", {
  set.seed(0); m <- 3; n <- 4
  x <- Variable(c(m, 1L), bounds = list(-2, 2))
  Y <- Variable(c(m, n), bounds = list(-1, 1))
  value(x) <- matrix(runif(m), m, 1)
  value(Y) <- matrix(runif(m * n), m, n)
  prob <- Problem(Minimize(sum_entries(x + Y)))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_true(all(abs(as.numeric(value(x)) - (-2)) < 1e-4))
  expect_true(all(abs(as.numeric(value(Y)) - (-1)) < 1e-4))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_affine_vector_atoms.py::TestAffineDiffEngine::test_index_stress
test_that("DE affine-vector: index stress sum(rows/cols/entries)", {
  set.seed(0); m <- 3; n <- 4
  X <- Variable(c(m, n), bounds = list(-2, 2))
  value(X) <- matrix(runif(m * n), m, n)
  ## CVXPY 0-based -> R 1-based: X[0,:] -> X[1,], X[:,2] -> X[,3],
  ## X[0,1] -> X[1,2], X[2,2] -> X[3,3].
  expr <- sum_entries(X[1, ]) + sum_entries(X[1, ]) +
    sum_entries(X[2, ]) + sum_entries(X[, 3]) + X[1, 2] + X[3, 3]
  prob <- Problem(Minimize(expr))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), -34.0, tolerance = 1e-4)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_affine_vector_atoms.py::TestAffineDiffEngine::test_duplicate_indices
test_that("DE affine-vector: duplicate fancy indices", {
  set.seed(0); m <- 3; n <- 3
  X <- Variable(c(m, n), bounds = list(-2, 2))
  value(X) <- matrix(runif(m * n), m, n)
  ## CVXPY X[[0,0],[1,1]] -> element pairs (X[0,1], X[0,1]); R SpecialIndex
  ## X[cbind(c(1,1), c(2,2))] -> c(X[1,2], X[1,2]). Expr collapses to sum(X).
  expr <- sum_entries(X[cbind(c(1L, 1L), c(2L, 2L))]) - 2 * X[1, 2] + sum_entries(X)
  prob <- Problem(Minimize(expr))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_true(all(abs(as.numeric(value(X)) - (-2)) < 1e-4))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_affine_vector_atoms.py::TestAffineDiffEngine::test_promote_row
test_that("DE affine-vector: promote scalar + row vector", {
  set.seed(0); n <- 4
  x <- Variable(bounds = list(-3, 3))
  Y <- Variable(c(1L, n), bounds = list(-2, 2))
  value(x) <- 2.0
  value(Y) <- matrix(runif(n), 1, n)
  prob <- Problem(Minimize(sum_entries(x + Y)))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), -3, tolerance = 1e-4)
  expect_true(all(abs(as.numeric(value(Y)) - (-2)) < 1e-4))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_affine_vector_atoms.py::TestAffineDiffEngine::test_promote_col
test_that("DE affine-vector: promote scalar + column vector", {
  set.seed(0); m <- 4
  x <- Variable(bounds = list(-3, 3))
  Y <- Variable(c(m, 1L), bounds = list(-2, 2))
  value(x) <- 2.0
  value(Y) <- matrix(runif(m), m, 1)
  prob <- Problem(Minimize(sum_entries(x + Y)))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), -3, tolerance = 1e-4)
  expect_true(all(abs(as.numeric(value(Y)) - (-2)) < 1e-4))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_affine_vector_atoms.py::TestAffineDiffEngine::test_promote_add
test_that("DE affine-vector: promote scalar + matrix", {
  set.seed(0)
  x <- Variable(bounds = list(-1, 1))
  Y <- Variable(c(2L, 2L), bounds = list(0, 2))
  value(x) <- 0.0
  value(Y) <- matrix(runif(4), 2, 2)
  prob <- Problem(Minimize(sum_entries(x + Y)))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), -1, tolerance = 1e-4)
  expect_true(all(abs(as.numeric(value(Y)) - 0) < 1e-4))
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_affine_vector_atoms.py::TestAffineDiffEngine::test_reshape
test_that("DE affine-vector: sum_squares(reshape(x, (4,2), F) - A)", {
  set.seed(0)
  x <- Variable(8, bounds = list(-5, 5))
  A <- matrix(runif(8), 4, 2)
  value(x) <- seq(-2, 2, length.out = 8)
  prob <- Problem(Minimize(sum_squares(reshape_expr(x, c(4L, 2L), order = "F") - A)))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), as.numeric(A), tolerance = 1e-4)  # column-major
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_affine_vector_atoms.py::TestAffineDiffEngine::test_broadcast
## R note: CVXPY subtracts a (8,1) column from a 1-D (8,) x, which numpy
## broadcasts to (8,8) (solution x_j = mean(A)). R's Variable(8) is already
## (8,1), so x - A is the (8,1) elementwise residual and the optimum is x = A.
## The diff engine is still exercised on the broadcast/reshape machinery.
test_that("DE affine-vector: sum_squares(x - A) broadcast", {
  set.seed(0)
  x <- Variable(8, bounds = list(-5, 5))
  A <- matrix(runif(8), 8, 1)
  value(x) <- seq(-2, 2, length.out = 8)
  prob <- Problem(Minimize(sum_squares(x - A)))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), as.numeric(A), tolerance = 1e-4)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_affine_vector_atoms.py::TestAffineDiffEngine::test_hstack
test_that("DE affine-vector: hstack of matmul residuals (NLP vs DCP)", {
  set.seed(0); m <- 5; n <- 3
  x <- Variable(c(n, 1L), bounds = list(-3, 3))
  y <- Variable(c(n, 1L), bounds = list(-2, 2))
  A1 <- matrix(runif(m * n), m, n); A2 <- matrix(runif(m * n), m, n)
  b1 <- matrix(runif(m), m, 1);     b2 <- matrix(runif(m), m, 1)
  blocks <- hstack(A1 %*% x + A2 %*% y - b1,
                   A1 %*% y + A2 %*% x - b2,
                   A2 %*% x - A1 %*% y)
  value(x) <- matrix(runif(n), n, 1); value(y) <- matrix(runif(n), n, 1)
  prob <- Problem(Minimize(sum_squares(blocks)))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  nlp_x <- value(x); nlp_y <- value(y); nlp_v <- value(prob)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), as.numeric(nlp_x), tolerance = 1e-4)
  expect_equal(as.numeric(value(y)), as.numeric(nlp_y), tolerance = 1e-4)
  expect_equal(value(prob), nlp_v, tolerance = 1e-4)
})

## @cvxpy nlp_tests/stress_tests_diff_engine/test_affine_vector_atoms.py::TestAffineDiffEngine::test_hstack_matrices
test_that("DE affine-vector: hstack of matrix matmul residuals (NLP vs DCP)", {
  set.seed(0); m <- 5; n <- 3
  X <- Variable(c(n, m), bounds = list(-3, 3))
  Y <- Variable(c(n, m), bounds = list(-2, 2))
  A1 <- matrix(runif(m * n), m, n); A2 <- matrix(runif(m * n), m, n)
  b1 <- matrix(runif(m * m), m, m); b2 <- matrix(runif(m * m), m, m)
  blocks <- hstack(A1 %*% X + A2 %*% Y - b1,
                   A1 %*% Y + A2 %*% X - b2,
                   A2 %*% X - A1 %*% Y)
  value(X) <- matrix(runif(n * m), n, m); value(Y) <- matrix(runif(n * m), n, m)
  prob <- Problem(Minimize(sum_squares(blocks)))
  nlp_derivative_check(prob)
  .de_nlp_solve(prob)
  nlp_x <- value(X); nlp_y <- value(Y); nlp_v <- value(prob)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(X)), as.numeric(nlp_x), tolerance = 1e-4)
  expect_equal(as.numeric(value(Y)), as.numeric(nlp_y), tolerance = 1e-4)
  expect_equal(value(prob), nlp_v, tolerance = 1e-4)
})
