## SCS Warm-Start Tests
## Tests for SCS warm-start support via solver cache (initial x/y/s).
## Mirrors CVXPY test_conic_solvers.py::test_warm_start (lines 293-304).

skip_if_not_installed("scs", minimum_version = "3.0")

# ── CVXPY parity: test_warm_start ────────────────────────────────────────────
## CVXPY test: solve exp-cone problem, then warm-start same problem.

## @cvxpy test_conic_solvers.py::TestSCS::test_warm_start
test_that("SCS warm-start matches CVXPY test_warm_start pattern", {
  x <- Variable(10)
  prob <- Problem(Minimize(sum_entries(exp(x))), list(sum_entries(x) == 1))

  r1 <- psolve(prob, solver = "SCS", warm_start = FALSE)
  expect_equal(status(prob), "optimal")

  r2 <- psolve(prob, solver = "SCS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r1, r2, tolerance = 1e-2)
})

# ── Basic warm-start ─────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("SCS warm-start: solve twice on same problem", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  r1 <- psolve(prob, solver = "SCS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r1, 2.0, tolerance = 1e-3)

  r2 <- psolve(prob, solver = "SCS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r2, 2.0, tolerance = 1e-3)
})

# ── Parameter change ─────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("SCS warm-start: parameter change between solves", {
  n <- 5L
  x <- Variable(n)
  b <- Parameter(n)
  prob <- Problem(Minimize(sum_squares(x - b)), list(x >= 0))

  value(b) <- rep(1, n)
  r1 <- psolve(prob, solver = "SCS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), rep(1, n), tolerance = 1e-2)

  value(b) <- c(2, 3, 4, 5, 6)
  r2 <- psolve(prob, solver = "SCS", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(2, 3, 4, 5, 6), tolerance = 1e-2)
})

# ── Cold vs warm correctness ────────────────────────────────────────────────

## @cvxpy NONE
test_that("SCS warm-start: cold and warm produce same solution", {
  m <- 50L
  n <- 20L
  set.seed(42L)
  A_mat <- matrix(rnorm(m * n), nrow = m, ncol = n)
  b_val <- rnorm(m)

  x <- Variable(n)
  prob <- Problem(Minimize(sum_squares(A_mat %*% x - b_val)))

  r_cold <- psolve(prob, solver = "SCS", warm_start = FALSE)
  r_warm <- psolve(prob, solver = "SCS", warm_start = TRUE)
  expect_equal(r_cold, r_warm, tolerance = 1e-2)
})

# ── Cold after warm ──────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("SCS cold solve after warm solve works", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  r1 <- psolve(prob, solver = "SCS", warm_start = TRUE)
  expect_equal(r1, 2.0, tolerance = 1e-3)

  r2 <- psolve(prob, solver = "SCS", warm_start = FALSE)
  expect_equal(r2, 2.0, tolerance = 1e-3)
})

# ── Independent caches per Problem ───────────────────────────────────────────

## @cvxpy NONE
test_that("SCS warm-start: two Problems have independent caches", {
  x1 <- Variable(2)
  prob1 <- Problem(Minimize(sum_entries(x1)), list(x1 >= 1))

  x2 <- Variable(2)
  prob2 <- Problem(Minimize(sum_entries(x2)), list(x2 >= 5))

  r1 <- psolve(prob1, solver = "SCS", warm_start = TRUE)
  r2 <- psolve(prob2, solver = "SCS", warm_start = TRUE)

  expect_equal(r1, 2.0, tolerance = 1e-3)
  expect_equal(r2, 10.0, tolerance = 1e-3)

  r1b <- psolve(prob1, solver = "SCS", warm_start = TRUE)
  expect_equal(r1b, 2.0, tolerance = 1e-3)
})
