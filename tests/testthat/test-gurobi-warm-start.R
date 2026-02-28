## Gurobi Warm-Start Tests
## Tests for Gurobi warm-start support via solver cache (model$start).
## Mirrors CVXPY test_conic_solvers.py::test_gurobi_warm_start and
## test_qp_solvers.py::test_gurobi_warmstart patterns.

skip_if_not_installed("gurobi")

# ── Basic warm-start (QP path) ──────────────────────────────────────────────

## @cvxpy NONE
test_that("Gurobi warm-start: solve twice on same QP", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x)), list(x >= 1))

  r1 <- psolve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r1, 2.0, tolerance = 1e-4)

  r2 <- psolve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(r2, 2.0, tolerance = 1e-4)
})

# ── Parameter change (QP path) ──────────────────────────────────────────────

## @cvxpy NONE
test_that("Gurobi warm-start: parameter change between QP solves", {
  n <- 5L
  x <- Variable(n)
  b <- Parameter(n)
  prob <- Problem(Minimize(sum_squares(x - b)), list(x >= 0))

  value(b) <- rep(1, n)
  r1 <- psolve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), rep(1, n), tolerance = 1e-3)

  value(b) <- c(2, 3, 4, 5, 6)
  r2 <- psolve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(2, 3, 4, 5, 6), tolerance = 1e-3)
})

# ── Cold vs warm correctness (QP path) ──────────────────────────────────────

## @cvxpy NONE
test_that("Gurobi warm-start: cold and warm produce same QP solution", {
  m <- 50L
  n <- 20L
  set.seed(42L)
  A_mat <- matrix(rnorm(m * n), nrow = m, ncol = n)
  b_param <- Parameter(m)

  x <- Variable(n)
  prob <- Problem(Minimize(sum_squares(A_mat %*% x - b_param)))

  value(b_param) <- rnorm(m)
  r_cold <- psolve(prob, solver = "GUROBI", warm_start = FALSE)
  r_warm <- psolve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(r_cold, r_warm, tolerance = 1e-3)

  ## Change parameter and test warm vs cold
  value(b_param) <- rnorm(m)
  r_warm2 <- psolve(prob, solver = "GUROBI", warm_start = TRUE)
  r_cold2 <- psolve(prob, solver = "GUROBI", warm_start = FALSE)
  expect_equal(r_warm2, r_cold2, tolerance = 1e-3)
})

# ── CVXPY parity: conic warm-start with parameter changes ───────────────────
## Mirrors test_conic_solvers.py::test_gurobi_warm_start (lines 1748-1809):
## parametric LP, change A/b, h, and c between warm-start solves.

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_warm_start
test_that("Gurobi warm-start: conic path with parameter changes (CVXPY parity)", {
  x <- Variable(2)
  A_par <- Parameter(c(2, 2))
  b_par <- Parameter(2)
  h_par <- Parameter(2)
  c_par <- Parameter(2)

  value(A_par) <- matrix(c(1, 0, 0, 0), 2, 2)
  value(b_par) <- c(1, 0)
  value(h_par) <- c(2, 2)
  value(c_par) <- c(1, 1)

  objective <- Maximize(c_par[1] * x[1] + c_par[2] * x[2])
  constraints <- list(x[1] <= h_par[1], x[2] <= h_par[2], A_par %*% x == b_par)
  prob <- Problem(objective, constraints)

  r1 <- psolve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(r1, 3, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 2), tolerance = 1e-3)

  ## Change A and b
  value(A_par) <- matrix(c(0, 0, 0, 1), 2, 2)
  value(b_par) <- c(0, 1)
  r2 <- psolve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(r2, 3, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(2, 1), tolerance = 1e-3)

  ## Change h
  value(A_par) <- matrix(c(1, 0, 0, 0), 2, 2)
  value(b_par) <- c(1, 0)
  value(h_par) <- c(1, 1)
  r3 <- psolve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(r3, 2, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-3)

  ## Change c
  value(h_par) <- c(2, 2)
  value(c_par) <- c(2, 1)
  r4 <- psolve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(r4, 4, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 2), tolerance = 1e-3)
})

# ── Cold after warm ──────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Gurobi cold solve after warm solve works", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  r1 <- psolve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(r1, 2.0, tolerance = 1e-4)

  r2 <- psolve(prob, solver = "GUROBI", warm_start = FALSE)
  expect_equal(r2, 2.0, tolerance = 1e-4)
})

# ── Independent caches per Problem ───────────────────────────────────────────

## @cvxpy NONE
test_that("Gurobi warm-start: two Problems have independent caches", {
  x1 <- Variable(2)
  prob1 <- Problem(Minimize(sum_entries(x1)), list(x1 >= 1))

  x2 <- Variable(2)
  prob2 <- Problem(Minimize(sum_entries(x2)), list(x2 >= 5))

  r1 <- psolve(prob1, solver = "GUROBI", warm_start = TRUE)
  r2 <- psolve(prob2, solver = "GUROBI", warm_start = TRUE)

  expect_equal(r1, 2.0, tolerance = 1e-4)
  expect_equal(r2, 10.0, tolerance = 1e-4)

  r1b <- psolve(prob1, solver = "GUROBI", warm_start = TRUE)
  expect_equal(r1b, 2.0, tolerance = 1e-4)
})
