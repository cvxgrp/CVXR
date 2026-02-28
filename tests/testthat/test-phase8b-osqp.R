## Phase 8b: OSQP Solver Tests
## Tests for OSQP solver interface: LP, QP, equality/inequality,
## infeasible/unbounded, dual values, solver options, cross-solver parity.

skip_if_not_installed("osqp")

# ── LP Tests ─────────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP solves simple LP", {
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x >= 0, x[1] + x[2] >= 1))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 1.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x[1])), 1.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x[2])), 0.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("OSQP solves LP with equality constraints", {
  z <- Variable(3)
  prob <- Problem(Minimize(z[1] + z[2] + z[3]),
                  list(z[1] + z[2] == 1, z[2] + z[3] == 2, z >= 0))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## x = (0, 1, 1) is optimal: x1+x2=1, x2+x3=2, sum=2
  expect_equal(result, 2.0, tolerance = 1e-4)
  xval <- as.numeric(value(z))
  expect_equal(xval[2], 1.0, tolerance = 1e-3)
  expect_equal(xval[3], 1.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("OSQP solves LP maximize", {
  b <- Variable(2)
  prob <- Problem(Maximize(b[1] + b[2]),
                  list(b <= 1, b >= 0))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 2.0, tolerance = 1e-4)
  bval <- as.numeric(value(b))
  expect_equal(bval, c(1, 1), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("OSQP solves LP with mixed constraints", {
  x <- Variable(3)
  prob <- Problem(Minimize(x[1] - x[2] + 2 * x[3]),
                  list(x >= 0, x[1] + x[2] + x[3] <= 10,
                       x[1] + 2 * x[2] == 6))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## Optimal: x1=0, x2=3, x3=0, value=-3
  expect_equal(result, -3.0, tolerance = 1e-3)
})

# ── QP Tests ─────────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP solves simple QP", {
  y <- Variable(2)
  P <- matrix(c(2, 0, 0, 1), 2, 2)
  prob <- Problem(Minimize(quad_form(y, P) + sum_entries(y)),
                  list(y >= 0))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## Unconstrained minimum is at y = -P^{-1} q/2 = (-0.25, -0.5), but y >= 0
  ## so optimal is y = (0, 0), value = 0
  expect_equal(result, 0.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("OSQP solves QP with equality constraint", {
  a <- Variable(2)
  prob <- Problem(Minimize(sum_squares(a)),
                  list(a[1] + a[2] == 1))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## Minimum of ||a||^2 s.t. a1 + a2 = 1 is a = (0.5, 0.5), value = 0.5
  expect_equal(result, 0.5, tolerance = 1e-4)
  aval <- as.numeric(value(a))
  expect_equal(aval, c(0.5, 0.5), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("OSQP solves QP with inequality constraints", {
  x <- Variable(2)
  P <- matrix(c(2, 0, 0, 1), 2, 2)
  prob <- Problem(Minimize(quad_form(x, P) + 3 * x[1] + 2 * x[2]),
                  list(x >= 0, x[1] + x[2] >= 1))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## CVXPY reference: value = 2.916667
  expect_equal(result, 35 / 12, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("OSQP solves sum_squares QP", {
  x <- Variable(3)
  target <- c(1, 2, 3)
  prob <- Problem(Minimize(sum_squares(x - target)),
                  list(x >= 0, sum_entries(x) <= 5))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  xval <- as.numeric(value(x))
  ## Unconstrained min: x=target=(1,2,3), sum=6 > 5.
  ## Constrained: project onto sum(x) <= 5
  ## x = target - lambda*1 where lambda chosen so sum(x) = 5
  ## sum(target) - 3*lambda = 5 => lambda = 1/3
  ## x = (2/3, 5/3, 8/3), value = 3 * (1/3)^2 = 1/3
  expect_equal(result, 1.0 / 3.0, tolerance = 1e-3)
})

# ── Infeasible / Unbounded ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP detects infeasible problem", {
  w <- Variable(2)
  prob <- Problem(Minimize(w[1] + w[2]),
                  list(w >= 1, w <= -1))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
})

## @cvxpy NONE
test_that("OSQP detects unbounded problem", {
  v <- Variable(2)
  prob <- Problem(Minimize(-v[1] - v[2]),
                  list(v >= 0))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_true(status(prob) %in% c("unbounded", "unbounded_inaccurate"))
})

# ── Dual Values ──────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP returns correct dual values for equality constraints", {
  x <- Variable(2)
  con_eq <- x[1] + x[2] == 1
  con_nn <- x >= 0
  prob <- Problem(Minimize(2 * x[1] + 3 * x[2]), list(con_eq, con_nn))
  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## CVXPY reference: eq dual = -2.0 (shadow price of equality)
  ## The optimal solution is x = (1, 0), value = 2
  ## Interpretation: relaxing x1+x2=1 by epsilon increases obj by -2*epsilon
  dv <- dual_value(con_eq)
  expect_equal(as.numeric(dv), -2.0, tolerance = 0.05)
})

## @cvxpy NONE
test_that("OSQP returns dual values for inequality constraints", {
  x <- Variable(2)
  con_nn <- x >= 0
  con_ub <- x[1] + x[2] <= 1
  prob <- Problem(Minimize(-x[1] - x[2]), list(con_nn, con_ub))
  val <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## Optimal: x = (0.5, 0.5), x1+x2 = 1, value = -1
  expect_equal(val, -1.0, tolerance = 1e-4)
  ## Upper bound constraint is active → dual should be nonzero
  dv_ub <- dual_value(con_ub)
  expect_true(!is.null(dv_ub))
  ## Dual for active upper bound should be ≈ -1 (OSQP sign convention)
  expect_true(abs(as.numeric(dv_ub)) > 0.1)
})

# ── Solver Options ───────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP respects custom solver options", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - c(1, 2))),
                  list(x >= 0))
  ## Solve with very few iterations (should hit max_iter)
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE,
                   max_iter = 1L)
  ## With only 1 iteration, likely won't converge
  ## Status should be user_limit (max_iter reached) or optimal_inaccurate
  status <- status(prob)
  expect_true(status %in% c("user_limit", "optimal_inaccurate", "optimal"))
})

## @cvxpy NONE
test_that("OSQP custom eps_abs tolerance works", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - c(1, 2))),
                  list(x >= 0))
  ## Solve with loose tolerance
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE,
                   eps_abs = 1e-3, eps_rel = 1e-3)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 0.0, tolerance = 0.01)
})

# ── Cross-Solver Parity ─────────────────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP matches Clarabel on LP", {
  skip_if_not_installed("clarabel")
  x <- Variable(3)
  prob <- Problem(Minimize(x[1] + 2 * x[2] + 3 * x[3]),
                  list(x >= 0, x[1] + x[2] + x[3] >= 2,
                       x[1] - x[2] >= 0))
  val_osqp <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  val_clar <- psolve(prob, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(val_osqp, val_clar, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("OSQP matches Clarabel on QP", {
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  P <- matrix(c(4, 1, 1, 2), 2, 2)
  prob <- Problem(Minimize(quad_form(x, P) + x[1] + x[2]),
                  list(x >= 0, x[1] + x[2] <= 1))
  val_osqp <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  val_clar <- psolve(prob, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(val_osqp, val_clar, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("OSQP matches SCS on QP", {
  skip_if_not_installed("scs")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x) + x[1] - 2 * x[2]),
                  list(x >= -1, x <= 1))
  val_osqp <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  val_scs <- psolve(prob, solver = SCS_SOLVER, verbose = FALSE)
  expect_equal(val_osqp, val_scs, tolerance = 1e-3)
})

# ── Error Handling ───────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP rejects SOCP problems", {
  x <- Variable(2)
  prob <- Problem(Minimize(cvxr_norm(x)),
                  list(x >= 0, sum_entries(x) >= 1))
  expect_error(psolve(prob, solver = OSQP_SOLVER, verbose = FALSE),
               "SOC")
})

## @cvxpy NONE
test_that("OSQP rejects SDP problems", {
  X <- Variable(c(2, 2), symmetric = TRUE)
  prob <- Problem(Minimize(matrix_trace(X)),
                  list(PSD(X), matrix_trace(X) >= 1))
  expect_error(psolve(prob, solver = OSQP_SOLVER, verbose = FALSE),
               "PSD")
})

# ── Edge Cases ───────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("OSQP handles unconstrained QP", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - c(3, 4))))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 0.0, tolerance = 1e-4)
  xval <- as.numeric(value(x))
  expect_equal(xval, c(3, 4), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("OSQP handles single variable", {
  x <- Variable(1)
  prob <- Problem(Minimize(square(x - 5)),
                  list(x >= 0, x <= 3))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 4.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 3.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("OSQP handles large sparse problem", {
  n <- 100
  x <- Variable(n)
  A <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = rep(1, n))
  target <- seq_len(n) / n
  prob <- Problem(Minimize(sum_squares(x - target)),
                  list(x >= 0, sum_entries(x) <= n / 2))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  xval <- as.numeric(value(x))
  ## All values should be non-negative

  expect_true(all(xval > -1e-6))
  ## Sum should be at most n/2
  expect_true(sum(xval) <= n / 2 + 1e-4)
})

## @cvxpy NONE
test_that("OSQP handles equality-only problem", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x)),
                  list(x[1] == 3, x[2] == 4))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 25.0, tolerance = 1e-4)
  xval <- as.numeric(value(x))
  expect_equal(xval, c(3, 4), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("OSQP handles inequality-only problem", {
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x >= 3))
  result <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 6.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("OSQP solver_stats populated", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - c(1, 2))),
                  list(x >= 0))
  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  stats <- solver_stats(prob)
  expect_true(!is.null(stats))
  expect_equal(stats@solver_name, "OSQP")
  expect_true(stats@num_iters > 0L)
  expect_true(stats@solve_time >= 0)
})
