## Regression tests for O(n^2) memory fix in quad canonicalizers.
## The quad canonicalizers (quad_over_lin_canon, power_canon) were
## calling as.matrix(Diagonal(n)) which densified an O(n) sparse
## identity into an O(n^2) dense matrix. Fixed by keeping P sparse
## as dgCMatrix throughout the pipeline.

## @cvxpy NONE -- R-specific memory fix, no CVXPY counterpart

test_that("OLS via sum_squares matches lm() at scale", {
  ## This is the original reporter's pattern (issues/marissa.R):
  ## sum_squares(y - X %*% beta) with large n_obs.
  set.seed(42)
  n <- 1000; k <- 20
  X <- matrix(rnorm(n * k), n, k)
  beta_true <- rnorm(k)
  y <- X %*% beta_true + rnorm(n, sd = 0.1)

  beta <- Variable(k)
  prob <- Problem(Minimize(sum_squares(y - X %*% beta)))
  psolve(prob)

  expect_equal(status(prob), "optimal")
  beta_ols <- unname(coef(lm(y ~ X + 0)))
  expect_equal(as.vector(value(beta)), beta_ols, tolerance = 1e-3)
})

test_that("power(x, 2) with large x produces correct result", {
  n <- 2000
  x <- Variable(n)
  prob <- Problem(Minimize(sum(power(x, 2))))
  psolve(prob)

  expect_equal(status(prob), "optimal")
  expect_equal(as.vector(value(x)), rep(0, n), tolerance = 1e-4)
})

test_that("quad_over_lin(large_expr, 1) produces correct result", {
  set.seed(42)
  n <- 2000; k <- 3
  A <- matrix(rnorm(n * k), n, k)
  x <- Variable(k)
  prob <- Problem(Minimize(quad_over_lin(A %*% x, 1)))
  psolve(prob)

  expect_equal(status(prob), "optimal")
  expect_equal(as.vector(value(x)), rep(0, k), tolerance = 1e-4)
})

test_that("sum_squares OLS works with both QP and conic solvers", {
  set.seed(42)
  n <- 500; k <- 5
  X <- matrix(rnorm(n * k), n, k)
  beta_true <- rnorm(k)
  y <- X %*% beta_true + rnorm(n, sd = 0.1)
  beta_ols <- unname(coef(lm(y ~ X + 0)))

  ## QP path (OSQP)
  beta1 <- Variable(k)
  prob1 <- Problem(Minimize(sum_squares(y - X %*% beta1)))
  psolve(prob1, solver = "OSQP")
  expect_equal(status(prob1), "optimal")
  expect_equal(as.vector(value(beta1)), beta_ols, tolerance = 1e-3)

  ## Conic path (Clarabel)
  beta2 <- Variable(k)
  prob2 <- Problem(Minimize(sum_squares(y - X %*% beta2)))
  psolve(prob2, solver = "CLARABEL")
  expect_equal(status(prob2), "optimal")
  expect_equal(as.vector(value(beta2)), beta_ols, tolerance = 1e-3)
})

test_that("quad_form with user-provided sparse P works correctly", {
  set.seed(42)
  n <- 100
  P <- Matrix::Diagonal(n, x = runif(n, 1, 5))
  ## Convert to dgCMatrix so coeff_quad_form takes sparse path
  P_sparse <- as(as(P, "generalMatrix"), "CsparseMatrix")
  x <- Variable(n)
  prob <- Problem(Minimize(quad_form(x, P_sparse)))
  psolve(prob)

  expect_equal(status(prob), "optimal")
  expect_equal(as.vector(value(x)), rep(0, n), tolerance = 1e-4)
})
