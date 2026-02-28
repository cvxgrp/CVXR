## OSQP Warm-Start Tests
## Tests for OSQP warm-start support via solver cache.
## Mirrors CVXPY test_qp_solvers.py::test_warm_start (lines 538-558).
## Requires osqp >= 1.0.0.

skip_if_not_installed("osqp", minimum_version = "1.0.0")

# ── Basic warm-start ─────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP warm-start: solve twice on same problem", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  r1 <- psolve(prob, solver = "OSQP", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r1, 2.0, tolerance = 1e-4)

  r2 <- psolve(prob, solver = "OSQP", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r2, 2.0, tolerance = 1e-4)
})

# ── Parameter change ─────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP warm-start: parameter change between solves", {
  n <- 5L
  x <- Variable(n)
  b <- Parameter(n)

  prob <- Problem(Minimize(sum_squares(x - b)), list(x >= 0))

  ## First solve: b = (1,1,1,1,1)
  value(b) <- rep(1, n)
  r1 <- psolve(prob, solver = "OSQP", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), rep(1, n), tolerance = 1e-3)

  ## Second solve: b = (2,3,4,5,6)
  value(b) <- c(2, 3, 4, 5, 6)
  r2 <- psolve(prob, solver = "OSQP", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(2, 3, 4, 5, 6), tolerance = 1e-3)
})

# ── Iteration improvement ───────────────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP warm-start: warm solve uses fewer or equal iterations", {
  n <- 20L
  x <- Variable(n)
  b <- Parameter(n)
  prob <- Problem(Minimize(sum_squares(x - b)), list(x >= 0))

  ## Cold solve
  value(b) <- seq_len(n)
  psolve(prob, solver = "OSQP", warm_start = TRUE)
  cold_iters <- solver_stats(prob)@num_iters

  ## Warm solve with tiny perturbation — should converge faster
  value(b) <- seq_len(n) + 0.001
  psolve(prob, solver = "OSQP", warm_start = TRUE)
  warm_iters <- solver_stats(prob)@num_iters

  expect_lte(warm_iters, cold_iters)
})

# ── Correctness: cold vs warm ────────────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP warm-start: cold and warm produce same solution", {
  x <- Variable(3)
  con <- list(x[1] + x[2] == 1, x[2] + x[3] == 2, x >= 0)

  ## Cold solve
  prob_cold <- Problem(Minimize(sum_squares(x)), con)
  r_cold <- psolve(prob_cold, solver = "OSQP", warm_start = FALSE)
  x_cold <- as.numeric(value(x))

  ## Warm solve (first solve builds cache, second uses it)
  x2 <- Variable(3)
  con2 <- list(x2[1] + x2[2] == 1, x2[2] + x2[3] == 2, x2 >= 0)
  prob_warm <- Problem(Minimize(sum_squares(x2)), con2)
  psolve(prob_warm, solver = "OSQP", warm_start = TRUE)
  r_warm <- psolve(prob_warm, solver = "OSQP", warm_start = TRUE)
  x_warm <- as.numeric(value(x2))

  expect_equal(r_cold, r_warm, tolerance = 1e-3)
  expect_equal(x_cold, x_warm, tolerance = 1e-3)
})

# ── Cold after warm ──────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP cold solve after warm solve creates fresh model", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  ## Warm solve (builds cache)
  r1 <- psolve(prob, solver = "OSQP", warm_start = TRUE)
  expect_equal(r1, 2.0, tolerance = 1e-4)

  ## Cold solve (should ignore cache)
  r2 <- psolve(prob, solver = "OSQP", warm_start = FALSE)
  expect_equal(r2, 2.0, tolerance = 1e-4)
})

# ── Independent caches per Problem ───────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP warm-start: two Problems have independent caches", {
  x1 <- Variable(2)
  prob1 <- Problem(Minimize(sum_entries(x1)), list(x1 >= 1))

  x2 <- Variable(2)
  prob2 <- Problem(Minimize(sum_entries(x2)), list(x2 >= 5))

  r1 <- psolve(prob1, solver = "OSQP", warm_start = TRUE)
  r2 <- psolve(prob2, solver = "OSQP", warm_start = TRUE)

  expect_equal(r1, 2.0, tolerance = 1e-4)
  expect_equal(r2, 10.0, tolerance = 1e-4)

  ## Solve prob1 again (warm) — should NOT be affected by prob2's cache
  r1b <- psolve(prob1, solver = "OSQP", warm_start = TRUE)
  expect_equal(r1b, 2.0, tolerance = 1e-4)
})

# ── Cross-solver unaffected ──────────────────────────────────────────────────

## @cvxpy NONE
test_that("warm_start=TRUE with non-OSQP solvers works (ignored gracefully)", {
  skip_if_not_installed("scs")
  skip_if_not_installed("clarabel")

  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  ## SCS with warm_start=TRUE (should be ignored, no error)
  r_scs <- psolve(prob, solver = "SCS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r_scs, 2.0, tolerance = 1e-3)

  ## Clarabel with warm_start=TRUE (should be ignored, no error)
  x2 <- Variable(2)
  prob2 <- Problem(Minimize(sum_entries(x2)), list(x2 >= 1))
  r_clar <- psolve(prob2, solver = "CLARABEL", warm_start = TRUE)
  expect_equal(status(prob2), "optimal")
  expect_equal(r_clar, 2.0, tolerance = 1e-3)
})

# ── QP with dual values ─────────────────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP warm-start preserves dual values", {
  x <- Variable(2)
  eq_con <- x[1] + x[2] == 1
  ineq_con <- x >= 0
  prob <- Problem(Minimize(sum_squares(x)), list(eq_con, ineq_con))

  ## First solve
  psolve(prob, solver = "OSQP", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(0.5, 0.5), tolerance = 1e-3)

  ## Second solve (warm)
  ## Use fresh Problem/Constraints to avoid cache leakage
  x2 <- Variable(2)
  eq_con2 <- x2[1] + x2[2] == 1
  ineq_con2 <- x2 >= 0
  prob2 <- Problem(Minimize(sum_squares(x2)), list(eq_con2, ineq_con2))

  psolve(prob2, solver = "OSQP", warm_start = TRUE)
  r2 <- psolve(prob2, solver = "OSQP", warm_start = TRUE)
  expect_equal(status(prob2), "optimal")
  expect_equal(as.numeric(value(x2)), c(0.5, 0.5), tolerance = 1e-3)
  ## Equality dual should be -1.0 (verified via CVXPY)
  expect_equal(as.numeric(dual_value(eq_con2)), -1.0, tolerance = 1e-2)
})

# ── CVXPY parity: test_warm_start ────────────────────────────────────────────
## Mirrors CVXPY test_qp_solvers.py::test_warm_start (lines 538-558):
## Parametric least-squares, cold-vs-warm and warm-vs-cold with parameter change.

## @cvxpy test_qp_solvers.py::TestQp::test_warm_start
test_that("OSQP warm-start matches CVXPY test_warm_start pattern", {
  m <- 200L
  n <- 100L
  set.seed(1L)
  A_mat <- matrix(rnorm(m * n), nrow = m, ncol = n)
  b <- Parameter(m)

  x <- Variable(n)
  prob <- Problem(Minimize(sum_squares(A_mat %*% x - b)))

  ## Cycle 1: same b — cold then warm should match
  value(b) <- rnorm(m)
  r1_cold <- psolve(prob, solver = "OSQP", warm_start = FALSE)
  r1_warm <- psolve(prob, solver = "OSQP", warm_start = TRUE)
  expect_equal(r1_cold, r1_warm, tolerance = 1e-3)

  ## Cycle 2: new b — warm then cold should match
  value(b) <- rnorm(m)
  r2_warm <- psolve(prob, solver = "OSQP", warm_start = TRUE)
  r2_cold <- psolve(prob, solver = "OSQP", warm_start = FALSE)
  expect_equal(r2_warm, r2_cold, tolerance = 1e-3)
})
