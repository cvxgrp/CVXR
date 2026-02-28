## Phase 8b: OSQP CVXPY Parity Tests
## Each test verified against CVXPY 1.8.1 with OSQP solver.
## Reference values obtained via: uv run python -c "import cvxpy as cp; ..."

skip_if_not_installed("osqp")

# ── quad_form_bound (CVXPY test_qp_solvers.py quad_form_bound) ───────────────

## @cvxpy test_qp_solvers.py::TestQp::test_all_solvers
test_that("OSQP parity: quad_form_bound", {
  ## CVXPY: P=[[13,12,-2],[12,17,6],[-2,6,12]], q=[[-22],[-14.5],[13]], r=1
  ## Minimize 0.5*quad_form(y, P) + t(q) %*% y + r, -1 <= y <= 1
  ## Expected: y = (1, 0.5, -1), value = -21.625
  P <- matrix(c(13, 12, -2, 12, 17, 6, -2, 6, 12), 3, 3)
  q <- c(-22, -14.5, 13)
  r <- 1
  y <- Variable(3)
  prob <- Problem(Minimize(0.5 * quad_form(y, P) + sum_entries(q * y) + r),
                  list(y >= -1, y <= 1))
  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), -21.625, tolerance = 1e-3)
  yval <- as.numeric(value(y))
  expect_equal(yval, c(1.0, 0.5, -1.0), tolerance = 1e-3)
})

# ── power (CVXPY test_qp_solvers.py power) ──────────────────────────────────

## @cvxpy test_qp_solvers.py::TestQp::test_all_solvers
test_that("OSQP parity: power/sum_squares unconstrained", {
  ## Minimize sum(power(x, 2)), unconstrained
  ## Expected: x = (0, 0), value = 0
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(power(x, 2))))
  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 0.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(0, 0), tolerance = 1e-3)
})

# ── huber_small (CVXPY test_qp_solvers.py huber_small) ──────────────────────

## @cvxpy test_qp_solvers.py::TestQp::test_all_solvers
test_that("OSQP parity: huber_small", {
  ## Minimize sum(huber(x)), x[3] >= 3
  ## Expected: x = (0, 0, 3), value = 5
  x <- Variable(3)
  prob <- Problem(Minimize(sum_entries(huber(x))),
                  list(x[3] >= 3))
  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 5.0, tolerance = 1e-3)
  xval <- as.numeric(value(x))
  expect_equal(xval[3], 3.0, tolerance = 1e-3)
  expect_equal(xval[1], 0.0, tolerance = 1e-2)
  expect_equal(xval[2], 0.0, tolerance = 1e-2)
})

# ── sparse_quad_form (CVXPY test_quad_form.py) ──────────────────────────────

## @cvxpy NONE
test_that("OSQP parity: sparse quad_form", {
  ## quad_form(x, I) with x == (1, 2)
  ## Expected: value = 5
  Q <- diag(2)
  x <- Variable(2)
  prob <- Problem(Minimize(quad_form(x, Q)),
                  list(x == c(1, 2)))
  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 5.0, tolerance = 1e-4)
})

# ── affine_problem (LP through OSQP, CVXPY test_qp_solvers.py) ─────────────

## @cvxpy test_qp_solvers.py::TestQp::test_all_solvers
test_that("OSQP parity: affine LP", {
  ## Minimize sum(x), x >= 0, A %*% x <= b (A, b positive)
  ## Expected: x = (0, 0), value = 0
  set.seed(0)
  A <- abs(matrix(rnorm(10), 5, 2))
  b <- abs(rnorm(5))
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)),
                  list(x >= 0, A %*% x <= b))
  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 0.0, tolerance = 1e-3)
  xval <- as.numeric(value(x))
  expect_equal(xval, c(0, 0), tolerance = 1e-3)
})

# ── maximize LP (CVXPY test_qp_solvers.py maximize_problem) ─────────────────

## @cvxpy test_qp_solvers.py::TestQp::test_all_solvers
test_that("OSQP parity: maximize LP", {
  ## Maximize -sum(x), x >= 0, A %*% x <= b (same A, b)
  ## Expected: x = (0, 0), value = 0
  set.seed(0)
  A <- abs(matrix(rnorm(10), 5, 2))
  b <- abs(rnorm(5))
  x <- Variable(2)
  prob <- Problem(Maximize(-sum_entries(x)),
                  list(x >= 0, A %*% x <= b))
  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 0.0, tolerance = 1e-3)
})

# ── sum_squares with equality (verified against CVXPY) ──────────────────────

## @cvxpy NONE
test_that("OSQP parity: sum_squares with equality", {
  ## Minimize sum_squares(x), x[1] + 2*x[2] == 3
  ## Expected: x = (0.6, 1.2), value = 1.8
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x)),
                  list(x[1] + 2 * x[2] == 3))
  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 1.8, tolerance = 1e-3)
  xval <- as.numeric(value(x))
  expect_equal(xval, c(0.6, 1.2), tolerance = 1e-3)
})

# ── quad_form with identity P (factor of 2 check) ──────────────────────────

## @cvxpy NONE
test_that("OSQP parity: quad_form identity + linear", {
  ## Minimize quad_form(x, I) + c(1,2,3) %*% x, x >= -5
  ## Expected: x = (-0.5, -1.0, -1.5), value = -3.5
  x <- Variable(3)
  c_vec <- c(1, 2, 3)
  prob <- Problem(Minimize(quad_form(x, diag(3)) + sum_entries(c_vec * x)),
                  list(x >= -5))
  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), -3.5, tolerance = 1e-3)
  xval <- as.numeric(value(x))
  expect_equal(xval, c(-0.5, -1.0, -1.5), tolerance = 1e-3)
})

# ── QP with mixed eq + ineq constraints ─────────────────────────────────────

## @cvxpy NONE
test_that("OSQP parity: QP with mixed eq and ineq", {
  ## Minimize sum_squares(x - (1,2,3)), sum(x) == 3, x >= 0
  ## Expected: x = (0, 1, 2), value = 3.0
  x <- Variable(3)
  target <- c(1, 2, 3)
  prob <- Problem(Minimize(sum_squares(x - target)),
                  list(x[1] + x[2] + x[3] == 3, x >= 0))
  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 3.0, tolerance = 1e-3)
  xval <- as.numeric(value(x))
  expect_equal(xval[1], 0.0, tolerance = 1e-2)
  expect_equal(xval[2], 1.0, tolerance = 1e-2)
  expect_equal(xval[3], 2.0, tolerance = 1e-2)
})

# ── square_affine / least-squares (CVXPY test_qp_solvers.py) ────────────────

## @cvxpy test_qp_solvers.py::TestQp::test_all_solvers
test_that("OSQP parity: least-squares", {
  ## Minimize sum_squares(A %*% x - b), unconstrained
  ## Expected: matches R's lstsq solution
  set.seed(42)
  m <- 10
  n <- 2
  A <- matrix(rnorm(m * n), m, n)
  b <- rnorm(m)
  x <- Variable(n)
  prob <- Problem(Minimize(sum_squares(A %*% x - b)))
  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## Compare with R's least-squares solution
  x_lstsq <- as.numeric(qr.solve(A, b))
  xval <- as.numeric(value(x))
  expect_equal(xval, x_lstsq, tolerance = 1e-3)
})

# ── Cross-solver parity: OSQP vs Clarabel on QP ─────────────────────────────

## @cvxpy NONE
test_that("OSQP parity: matches Clarabel on quad_form_bound", {
  skip_if_not_installed("clarabel")
  P <- matrix(c(13, 12, -2, 12, 17, 6, -2, 6, 12), 3, 3)
  q <- c(-22, -14.5, 13)
  y <- Variable(3)
  prob <- Problem(Minimize(0.5 * quad_form(y, P) + sum_entries(q * y) + 1),
                  list(y >= -1, y <= 1))
  val_osqp <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  y_osqp <- as.numeric(value(y))
  val_clar <- psolve(prob, solver = CLARABEL_SOLVER, verbose = FALSE)
  y_clar <- as.numeric(value(y))
  expect_equal(val_osqp, val_clar, tolerance = 1e-3)
  expect_equal(y_osqp, y_clar, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("OSQP parity: matches SCS on huber", {
  skip_if_not_installed("scs")
  x <- Variable(3)
  prob <- Problem(Minimize(sum_entries(huber(x))),
                  list(x[3] >= 3))
  val_osqp <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  val_scs <- psolve(prob, solver = SCS_SOLVER, verbose = FALSE)
  expect_equal(val_osqp, val_scs, tolerance = 1e-2)
})

# ── Error handling parity ───────────────────────────────────────────────────

## @cvxpy test_qp_solvers.py::TestQpSolverValidation::test_qp_solver_rejects_soc_cones
test_that("OSQP parity: rejects SOCP problem", {
  ## CVXPY: test_invalid_solver — lambda_max needs SDP, not QP
  ## In CVXR: norm() creates SOC constraint, OSQP should reject
  x <- Variable(3)
  prob <- Problem(Minimize(cvxr_norm(x, 2)),
                  list(sum_entries(x) == 1))
  expect_error(psolve(prob, solver = OSQP_SOLVER, verbose = FALSE),
               "SOC")
})

# ── solver_stats (CVXPY test_problem.py test_solver_stats) ──────────────────

## @cvxpy test_problem.py::TestProblem::test_solver_stats
test_that("OSQP parity: solver_stats populated", {
  ## CVXPY: stats.solve_time > 0, stats.num_iters > 0
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)),
                  list(x == 0))
  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  stats <- solver_stats(prob)
  expect_equal(stats@solver_name, "OSQP")
  expect_true(stats@solve_time >= 0)
  expect_true(stats@num_iters > 0L)
})

# ── Equivalent forms: sum_squares vs quad_form ──────────────────────────────

## @cvxpy test_qp_solvers.py::TestQp::test_all_solvers
test_that("OSQP parity: sum_squares == quad_form equivalent", {
  ## CVXPY test_qp_solvers.py: equivalent_forms 1 vs 2
  ## sum_squares(A %*% x - b) = quad_form(x, A'A) - 2*b'Ax + b'b
  set.seed(42)
  m <- 10
  n <- 2
  A <- matrix(rnorm(m * n), m, n)
  b <- rnorm(m)
  x1 <- Variable(n)
  prob1 <- Problem(Minimize(sum_squares(A %*% x1 - b)),
                   list(x1[1] + x1[2] >= 0))
  val1 <- psolve(prob1, solver = OSQP_SOLVER, verbose = FALSE)

  P <- t(A) %*% A
  q <- -2 * t(A) %*% b
  r <- sum(b^2)
  x2 <- Variable(n)
  prob2 <- Problem(Minimize(quad_form(x2, P) + sum_entries(as.numeric(q) * x2) + r),
                   list(x2[1] + x2[2] >= 0))
  val2 <- psolve(prob2, solver = OSQP_SOLVER, verbose = FALSE)

  expect_equal(val1, val2, tolerance = 1e-3)
  expect_equal(as.numeric(value(x1)), as.numeric(value(x2)), tolerance = 1e-3)
})
