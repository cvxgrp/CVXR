## Phase 10a: HiGHS Solver Tests
## Tests for HiGHS solver interface: LP, QP, equality/inequality,
## infeasible/unbounded, dual values, solver options, cross-solver parity,
## error handling for unsupported cones, solver stats.

skip_if_not_installed("highs")

# ── LP Tests ─────────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("HiGHS solves simple LP", {
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x >= 0, x[1] + x[2] >= 1))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 1.0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x[1])), 1.0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x[2])), 0.0, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("HiGHS solves LP with equality constraints", {
  z <- Variable(3)
  prob <- Problem(Minimize(z[1] + z[2] + z[3]),
                  list(z[1] + z[2] == 1, z[2] + z[3] == 2, z >= 0))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## x = (0, 1, 1) is optimal: x1+x2=1, x2+x3=2, sum=2
  expect_equal(result, 2.0, tolerance = 1e-5)
  xval <- as.numeric(value(z))
  expect_equal(xval[2], 1.0, tolerance = 1e-5)
  expect_equal(xval[3], 1.0, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("HiGHS solves LP maximize", {
  b <- Variable(2)
  prob <- Problem(Maximize(b[1] + b[2]),
                  list(b <= 1, b >= 0))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 2.0, tolerance = 1e-5)
  bval <- as.numeric(value(b))
  expect_equal(bval, c(1, 1), tolerance = 1e-5)
})

## @cvxpy NONE
test_that("HiGHS solves LP with mixed constraints", {
  x <- Variable(3)
  prob <- Problem(Minimize(x[1] - x[2] + 2 * x[3]),
                  list(x >= 0, x[1] + x[2] + x[3] <= 10,
                       x[1] + 2 * x[2] == 6))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## Optimal: x1=0, x2=3, x3=0, value=-3
  expect_equal(result, -3.0, tolerance = 1e-5)
})

# ── QP Tests ─────────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("HiGHS solves simple QP", {
  y <- Variable(2)
  P <- matrix(c(2, 0, 0, 1), 2, 2)
  prob <- Problem(Minimize(quad_form(y, P) + sum_entries(y)),
                  list(y >= 0))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## Unconstrained minimum at y = -P^{-1}*1/2 = (-0.25, -0.5), but y >= 0
  ## so optimal is y = (0, 0), value = 0
  expect_equal(result, 0.0, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("HiGHS solves QP with equality constraint", {
  a <- Variable(2)
  prob <- Problem(Minimize(sum_squares(a)),
                  list(a[1] + a[2] == 1))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## Minimum of ||a||^2 s.t. a1 + a2 = 1 is a = (0.5, 0.5), value = 0.5
  expect_equal(result, 0.5, tolerance = 1e-5)
  aval <- as.numeric(value(a))
  expect_equal(aval, c(0.5, 0.5), tolerance = 1e-5)
})

## @cvxpy NONE
test_that("HiGHS solves QP with inequality constraints", {
  x <- Variable(2)
  P <- matrix(c(2, 0, 0, 1), 2, 2)
  prob <- Problem(Minimize(quad_form(x, P) + 3 * x[1] + 2 * x[2]),
                  list(x >= 0, x[1] + x[2] >= 1))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## CVXPY reference: value = 35/12 ≈ 2.916667
  expect_equal(result, 35 / 12, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("HiGHS solves sum_squares QP", {
  x <- Variable(3)
  target <- c(1, 2, 3)
  prob <- Problem(Minimize(sum_squares(x - target)),
                  list(x >= 0, sum_entries(x) <= 5))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## x = target - lambda*1 where lambda = 1/3, value = 1/3
  expect_equal(result, 1.0 / 3.0, tolerance = 1e-4)
})

# ── Infeasible / Unbounded ───────────────────────────────────────────────────

## @cvxpy NONE
test_that("HiGHS detects infeasible problem", {
  w <- Variable(2)
  prob <- Problem(Minimize(w[1] + w[2]),
                  list(w >= 1, w <= -1))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_true(status(prob) %in%
    c("infeasible", "infeasible_inaccurate", "infeasible_or_unbounded"))
})

## @cvxpy NONE
test_that("HiGHS detects unbounded problem", {
  v <- Variable(2)
  prob <- Problem(Minimize(-v[1] - v[2]),
                  list(v >= 0))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_true(status(prob) %in%
    c("unbounded", "unbounded_inaccurate", "infeasible_or_unbounded"))
})

# ── Dual Values ──────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("HiGHS returns correct dual values for equality constraints", {
  x <- Variable(2)
  con_eq <- x[1] + x[2] == 1
  con_nn <- x >= 0
  prob <- Problem(Minimize(2 * x[1] + 3 * x[2]), list(con_eq, con_nn))
  psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## Optimal: x = (1, 0), value = 2. Eq dual = -2.0
  dv <- dual_value(con_eq)
  expect_equal(as.numeric(dv), -2.0, tolerance = 0.05)
})

## @cvxpy NONE
test_that("HiGHS returns dual values for inequality constraints", {
  x <- Variable(2)
  con_nn <- x >= 0
  con_ub <- x[1] + x[2] <= 1
  prob <- Problem(Minimize(-x[1] - x[2]), list(con_nn, con_ub))
  val <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(val, -1.0, tolerance = 1e-5)
  dv_ub <- dual_value(con_ub)
  expect_true(!is.null(dv_ub))
  expect_true(abs(as.numeric(dv_ub)) > 0.1)
})

## @cvxpy NONE
test_that("HiGHS dual values match OSQP on LP", {
  skip_if_not_installed("osqp")
  x <- Variable(2)
  con_eq <- x[1] + x[2] == 2
  con_nn <- x >= 0
  prob <- Problem(Minimize(x[1] + 3 * x[2]), list(con_eq, con_nn))

  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  dv_osqp <- as.numeric(dual_value(con_eq))

  psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  dv_highs <- as.numeric(dual_value(con_eq))

  expect_equal(dv_highs, dv_osqp, tolerance = 0.05)
})

## @cvxpy NONE
test_that("HiGHS dual values match OSQP on QP", {
  skip_if_not_installed("osqp")
  x <- Variable(2)
  con_eq <- x[1] + x[2] == 1
  prob <- Problem(Minimize(sum_squares(x)), list(con_eq))

  psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  dv_osqp <- as.numeric(dual_value(con_eq))

  psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  dv_highs <- as.numeric(dual_value(con_eq))

  expect_equal(dv_highs, dv_osqp, tolerance = 0.05)
})

# ── Cross-Solver Parity ─────────────────────────────────────────────────────

## @cvxpy NONE
test_that("HiGHS matches Clarabel on LP", {
  skip_if_not_installed("clarabel")
  x <- Variable(3)
  prob <- Problem(Minimize(x[1] + 2 * x[2] + 3 * x[3]),
                  list(x >= 0, x[1] + x[2] + x[3] >= 2,
                       x[1] - x[2] >= 0))
  val_highs <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  val_clar <- psolve(prob, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(val_highs, val_clar, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("HiGHS matches Clarabel on QP", {
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  P <- matrix(c(4, 1, 1, 2), 2, 2)
  prob <- Problem(Minimize(quad_form(x, P) + x[1] + x[2]),
                  list(x >= 0, x[1] + x[2] <= 1))
  val_highs <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  val_clar <- psolve(prob, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(val_highs, val_clar, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("HiGHS matches OSQP on LP", {
  skip_if_not_installed("osqp")
  x <- Variable(3)
  prob <- Problem(Minimize(x[1] + 2 * x[2] + 3 * x[3]),
                  list(x >= 0, x[1] + x[2] + x[3] >= 2,
                       x[1] - x[2] >= 0))
  val_highs <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  val_osqp <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(val_highs, val_osqp, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("HiGHS matches OSQP on QP", {
  skip_if_not_installed("osqp")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x) + x[1] - 2 * x[2]),
                  list(x >= -1, x <= 1))
  val_highs <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  val_osqp <- psolve(prob, solver = OSQP_SOLVER, verbose = FALSE)
  expect_equal(val_highs, val_osqp, tolerance = 1e-3)
})

# ── Error Handling ───────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("HiGHS rejects SOCP problems", {
  x <- Variable(2)
  prob <- Problem(Minimize(cvxr_norm(x)),
                  list(x >= 0, sum_entries(x) >= 1))
  expect_error(psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE),
               "SOC")
})

## @cvxpy NONE
test_that("HiGHS rejects SDP problems", {
  X <- Variable(c(2, 2), symmetric = TRUE)
  prob <- Problem(Minimize(matrix_trace(X)),
                  list(PSD(X), matrix_trace(X) >= 1))
  expect_error(psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE),
               "PSD")
})

## @cvxpy NONE
test_that("HiGHS rejects ExpCone problems", {
  x <- Variable()
  y <- Variable()
  prob <- Problem(Minimize(x),
                  list(ExpCone(x, Constant(1), y), y <= 1, x >= 0))
  expect_error(psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE),
               "ExpCone")
})

# ── Edge Cases ───────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("HiGHS handles unconstrained QP", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - c(3, 4))))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 0.0, tolerance = 1e-5)
  xval <- as.numeric(value(x))
  expect_equal(xval, c(3, 4), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("HiGHS handles single variable", {
  x <- Variable(1)
  prob <- Problem(Minimize(square(x - 5)),
                  list(x >= 0, x <= 3))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 4.0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), 3.0, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("HiGHS handles large sparse LP", {
  n <- 100
  x <- Variable(n)
  target <- seq_len(n) / n
  prob <- Problem(Minimize(sum_entries(x)),
                  list(x >= target, sum_entries(x) <= 200))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  xval <- as.numeric(value(x))
  ## Optimal: x_i = target_i = i/n (all lower bounds active, sum=50.5 < 200)
  expect_equal(xval, target, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("HiGHS handles equality-only problem", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x)),
                  list(x[1] == 3, x[2] == 4))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 25.0, tolerance = 1e-5)
  xval <- as.numeric(value(x))
  expect_equal(xval, c(3, 4), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("HiGHS handles inequality-only problem", {
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(x >= 3))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 6.0, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("HiGHS solver_stats populated", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - c(1, 2))),
                  list(x >= 0))
  psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  stats <- solver_stats(prob)
  expect_true(!is.null(stats))
  expect_equal(stats@solver_name, "HIGHS")
})

# ── Solver constant exported ────────────────────────────────────────────────

## @cvxpy NONE
test_that("HIGHS_SOLVER constant is exported and correct", {
  expect_equal(HIGHS_SOLVER, "HIGHS")
  expect_true(exists("HIGHS_SOLVER"))
})

# ── Warm start ───────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("HiGHS warm start does not crash", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - c(1, 2))),
                  list(x >= 0))
  ## First solve
  psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## Second solve (warm_start=TRUE)
  psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE, warm_start = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(1, 2), tolerance = 1e-5)
})

# ── Nonneg attribute ────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("HiGHS handles nonneg attribute", {
  x <- Variable(2, nonneg = TRUE)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x[1] + x[2] >= 1))
  result <- psolve(prob, solver = HIGHS_SOLVER, verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 1.0, tolerance = 1e-5)
})
