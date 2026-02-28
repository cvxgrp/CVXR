## Clarabel Warm-Start Tests
## Tests for Clarabel warm-start support via persistent solver + solver_update().

skip_if_not_installed("clarabel", minimum_version = "0.11.2")

# ── Basic warm-start ─────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Clarabel warm-start: solve twice on same problem", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  r1 <- psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r1, 2.0, tolerance = 1e-4)

  r2 <- psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r2, 2.0, tolerance = 1e-4)
})

# ── Parameter change ─────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Clarabel warm-start: parameter change between solves", {
  n <- 5L
  x <- Variable(n)
  b <- Parameter(n)

  prob <- Problem(Minimize(sum_squares(x - b)), list(x >= 0))

  ## First solve: b = (1,1,1,1,1)
  value(b) <- rep(1, n)
  r1 <- psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), rep(1, n), tolerance = 1e-3)

  ## Second solve: b = (2,3,4,5,6)
  value(b) <- c(2, 3, 4, 5, 6)
  r2 <- psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(2, 3, 4, 5, 6), tolerance = 1e-3)
})

# ── Correctness: cold vs warm ────────────────────────────────────────────────

## @cvxpy NONE
test_that("Clarabel warm-start: cold and warm produce same solution", {
  x <- Variable(3)
  con <- list(x[1] + x[2] == 1, x[2] + x[3] == 2, x >= 0)

  ## Cold solve
  prob_cold <- Problem(Minimize(sum_squares(x)), con)
  r_cold <- psolve(prob_cold, solver = "CLARABEL", warm_start = FALSE)
  x_cold <- as.numeric(value(x))

  ## Warm solve (first solve builds cache, second uses it)
  x2 <- Variable(3)
  con2 <- list(x2[1] + x2[2] == 1, x2[2] + x2[3] == 2, x2 >= 0)
  prob_warm <- Problem(Minimize(sum_squares(x2)), con2)
  psolve(prob_warm, solver = "CLARABEL", warm_start = TRUE)
  r_warm <- psolve(prob_warm, solver = "CLARABEL", warm_start = TRUE)
  x_warm <- as.numeric(value(x2))

  expect_equal(r_cold, r_warm, tolerance = 1e-4)
  expect_equal(x_cold, x_warm, tolerance = 1e-4)
})

# ── Cold after warm ──────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Clarabel cold solve after warm solve creates fresh solver", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  ## Warm solve (builds cache)
  r1 <- psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  expect_equal(r1, 2.0, tolerance = 1e-4)

  ## Cold solve (should ignore cache)
  r2 <- psolve(prob, solver = "CLARABEL", warm_start = FALSE)
  expect_equal(r2, 2.0, tolerance = 1e-4)
})

# ── Independent caches per Problem ───────────────────────────────────────────

## @cvxpy NONE
test_that("Clarabel warm-start: two Problems have independent caches", {
  x1 <- Variable(2)
  prob1 <- Problem(Minimize(sum_entries(x1)), list(x1 >= 1))

  x2 <- Variable(2)
  prob2 <- Problem(Minimize(sum_entries(x2)), list(x2 >= 5))

  r1 <- psolve(prob1, solver = "CLARABEL", warm_start = TRUE)
  r2 <- psolve(prob2, solver = "CLARABEL", warm_start = TRUE)

  expect_equal(r1, 2.0, tolerance = 1e-4)
  expect_equal(r2, 10.0, tolerance = 1e-4)

  ## Solve prob1 again (warm) — should NOT be affected by prob2's cache
  r1b <- psolve(prob1, solver = "CLARABEL", warm_start = TRUE)
  expect_equal(r1b, 2.0, tolerance = 1e-4)
})

# ── QP with dual values ─────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Clarabel warm-start preserves dual values", {
  x <- Variable(2)
  eq_con <- x[1] + x[2] == 1
  ineq_con <- x >= 0
  prob <- Problem(Minimize(sum_squares(x)), list(eq_con, ineq_con))

  ## First solve
  psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(0.5, 0.5), tolerance = 1e-3)

  ## Second solve (warm) — fresh Problem to avoid cache leakage
  x2 <- Variable(2)
  eq_con2 <- x2[1] + x2[2] == 1
  ineq_con2 <- x2 >= 0
  prob2 <- Problem(Minimize(sum_squares(x2)), list(eq_con2, ineq_con2))

  psolve(prob2, solver = "CLARABEL", warm_start = TRUE)
  r2 <- psolve(prob2, solver = "CLARABEL", warm_start = TRUE)
  expect_equal(status(prob2), "optimal")
  expect_equal(as.numeric(value(x2)), c(0.5, 0.5), tolerance = 1e-3)
  ## Equality dual should be -1.0 (verified via CVXPY)
  expect_equal(as.numeric(dual_value(eq_con2)), -1.0, tolerance = 1e-2)
})

# ── CVXPY parity: test_clarabel_parameter_update ─────────────────────────────
## Mirrors CVXPY test_conic_solvers.py::test_clarabel_parameter_update
## (lines 442-475). Exercises all 4 parameter types (P, A, b, q) with
## equality + inequality constraints. Three cycles:
##   1. cold → warm (same params)
##   2. warm → cold (new params)
##   3. consecutive cold (no data update)

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_parameter_update
test_that("Clarabel warm-start: CVXPY parameter update parity", {
  set.seed(42L)

  x <- Variable(2)
  P_par <- Parameter(nonneg = TRUE)
  A_par <- Parameter(4)
  b_par <- Parameter(2, nonneg = TRUE)
  q_par <- Parameter(2)

  update_parameters <- function() {
    value(P_par) <<- runif(1)
    value(A_par) <<- rnorm(4)
    value(b_par) <<- runif(2)
    value(q_par) <<- rnorm(2)
  }

  Q_mat <- matrix(1, 2, 2)
  prob <- Problem(
    Minimize(P_par * power(x[1], 2) + quad_form(x, Q_mat) + t(q_par) %*% x),
    list(
      A_par[1] * x[1] + A_par[2] * x[2] == b_par[1],
      A_par[3] * x[1] + A_par[4] * x[2] <= b_par[2]
    )
  )

  ## Cycle 1: cold then warm — same parameters
  update_parameters()
  r1 <- psolve(prob, solver = "CLARABEL", warm_start = FALSE)
  r2 <- psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  expect_equal(r1, r2, tolerance = 1e-4)

  ## Cycle 2: warm then cold — new parameters
  update_parameters()
  r3 <- psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  r4 <- psolve(prob, solver = "CLARABEL", warm_start = FALSE)
  expect_equal(r3, r4, tolerance = 1e-4)

  ## Cycle 3: consecutive cold solve, no data update
  r5 <- psolve(prob, solver = "CLARABEL", warm_start = FALSE)
  expect_equal(r5, r4, tolerance = 1e-4)
})

# ── CVXPY parity: parametric least-squares ───────────────────────────────────
## Mirrors CVXPY test_qp_solvers.py::test_warm_start (lines 538-558),
## adapted for Clarabel (conic path).

## @cvxpy test_qp_solvers.py::TestQp::test_warm_start
test_that("Clarabel warm-start: parametric least-squares parity", {
  m <- 200L
  n <- 100L
  set.seed(1L)
  A_mat <- matrix(rnorm(m * n), nrow = m, ncol = n)
  b <- Parameter(m)

  x <- Variable(n)
  prob <- Problem(Minimize(sum_squares(A_mat %*% x - b)))

  ## Cycle 1: same b — cold then warm should match
  value(b) <- rnorm(m)
  r1_cold <- psolve(prob, solver = "CLARABEL", warm_start = FALSE)
  r1_warm <- psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  expect_equal(r1_cold, r1_warm, tolerance = 1e-3)

  ## Cycle 2: new b — warm then cold should match
  value(b) <- rnorm(m)
  r2_warm <- psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  r2_cold <- psolve(prob, solver = "CLARABEL", warm_start = FALSE)
  expect_equal(r2_warm, r2_cold, tolerance = 1e-3)
})

# ── Dimension change fallback ────────────────────────────────────────────────

## @cvxpy NONE
test_that("Clarabel warm-start falls back on dimension change", {
  x <- Variable(2)
  prob1 <- Problem(Minimize(sum_squares(x)), list(x >= 0))
  psolve(prob1, solver = "CLARABEL", warm_start = TRUE)

  ## Different problem structure — should fall back to cold
  y <- Variable(3)
  prob2 <- Problem(Minimize(sum_squares(y)), list(y >= 0))
  val <- psolve(prob2, solver = "CLARABEL", warm_start = TRUE)
  expect_true(is.finite(val))
  expect_equal(val, 0.0, tolerance = 1e-4)
})
