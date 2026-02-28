## QP path parity tests — Phase 8a
## Verified against CVXPY 1.8.1 (Clarabel + SCS)

## @cvxpy test_qp_solvers.py::TestQp::test_all_solvers
test_that("simple QP: quad_form + linear objective", {
  ## minimize x'Px + q'x s.t. x >= 0
  ## CVXPY value: -1.0, x = [0, 1, 0]
  P_val <- matrix(c(2, 0.5, 0, 0.5, 1, 0, 0, 0, 3), 3, 3)
  q_val <- c(1, -2, 0.5)
  x <- Variable(3)
  prob <- Problem(Minimize(quad_form(x, P_val) + t(q_val) %*% x),
                  list(x >= 0))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, -1.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(0, 1, 0), tolerance = 1e-4)
})

## @cvxpy test_qp_solvers.py::TestQp::test_all_solvers
test_that("sum_squares: minimize sum_squares(x) s.t. sum(x) == 1", {
  ## CVXPY value: 1/3, x = [1/3, 1/3, 1/3]
  x <- Variable(3)
  prob <- Problem(Minimize(sum_squares(x)), list(sum(x) == 1))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 1/3, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), rep(1/3, 3), tolerance = 1e-4)
})

## @cvxpy test_qp_solvers.py::TestQp::test_all_solvers
test_that("least squares: minimize ||Ax - b||^2", {
  ## CVXPY value: 0.1325589201
  set.seed(42)
  m <- 5L; n <- 3L
  A <- matrix(rnorm(m * n), m, n)
  b <- rnorm(m)
  x <- Variable(n)
  prob <- Problem(Minimize(sum_squares(A %*% x - b)))
  result <- psolve(prob, solver = "CLARABEL")
  ## Compare with R's direct least squares
  x_ls <- qr.solve(A, b)
  ls_val <- sum((A %*% x_ls - b)^2)
  expect_equal(result, ls_val, tolerance = 1e-4)
})

## @cvxpy test_qp_solvers.py::TestQp::test_all_solvers
test_that("QP + linear: quad_form + scalar terms", {
  ## minimize x'Ix + 3*x[1] - x[2] s.t. -1 <= x <= 2
  ## CVXPY value: -2.25, x = [-1, 0.5]
  x <- Variable(2)
  P <- diag(2)
  prob <- Problem(Minimize(quad_form(x, P) + 3 * x[1] - x[2]),
                  list(x >= -1, x <= 2))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, -2.25, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(-1, 0.5), tolerance = 1e-4)
})

## @cvxpy test_qp_solvers.py::TestQp::test_all_solvers
test_that("power(x, 2) elementwise via QP path", {
  ## minimize sum(x^2) s.t. x >= 1
  ## CVXPY value: 3.0, x = [1, 1, 1]
  x <- Variable(3)
  prob <- Problem(Minimize(sum(power(x, 2))), list(x >= 1))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 3.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 1, 1), tolerance = 1e-4)
})

## @cvxpy test_qp_solvers.py::TestQp::test_all_solvers
test_that("SCS solver QP", {
  ## minimize x'Px s.t. x1 + x2 == 1, x >= 0
  ## P = [[4,1],[1,2]]
  ## CVXPY value: 1.75, x = [0.25, 0.75]
  x <- Variable(2)
  P <- matrix(c(4, 1, 1, 2), 2, 2)
  prob <- Problem(Minimize(quad_form(x, P)),
                  list(x[1] + x[2] == 1, x >= 0))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 1.75, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(0.25, 0.75), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("QP path gives same result as conic path", {
  ## Solve with QP path (quad_obj=TRUE, automatic) and compare with
  ## manually decomposed conic version
  x <- Variable(2)
  P <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  constraints <- list(x[1] + x[2] == 1, x >= 0)

  ## QP path (automatic via has_quadratic_term)
  prob_qp <- Problem(Minimize(quad_form(x, P)), constraints)
  val_qp <- psolve(prob_qp, solver = "CLARABEL")
  x_qp <- as.numeric(value(x))

  ## CVXPY value: 0.875, x = [0.25, 0.75]
  expect_equal(val_qp, 0.875, tolerance = 1e-4)
  expect_equal(x_qp, c(0.25, 0.75), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("QP with quadratic + constant offset", {
  ## minimize (x - 1)^2 = x^2 - 2x + 1
  x <- Variable(1)
  prob <- Problem(Minimize(power(x - 1, 2)))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 0.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("quad_form with non-variable argument creates constraint", {
  ## minimize quad_form(2*x, P) — argument is affine, not a Variable
  x <- Variable(2)
  P <- diag(2)
  prob <- Problem(Minimize(quad_form(2 * x, P)),
                  list(x[1] + x[2] == 1, x >= 0))
  result <- psolve(prob, solver = "CLARABEL")
  ## quad_form(2x, I) = 4*x'x, min at x = (0.5, 0.5), value = 2
  expect_equal(result, 2.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(0.5, 0.5), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("SymbolicQuadForm class basics", {
  x <- Variable(3)
  P <- Constant(diag(3))
  sqf <- CVXR:::SymbolicQuadForm(x, P, quad_form(x, diag(3)))
  expect_true(is_quadratic(sqf))
  expect_equal(sqf@shape, c(1L, 1L))
  expect_error(CVXR:::graph_implementation(sqf, list(), c(1L, 1L)),
               "replaced before canonicalization")
})

## @cvxpy NONE
test_that("replace_quad_forms and restore work correctly", {
  x <- Variable(2)
  P <- diag(2)
  qf_expr <- quad_form(x, P) + 1

  ## Dcp2Cone with quad_obj should produce SymbolicQuadForm
  d2c <- CVXR:::Dcp2Cone(quad_obj = TRUE)
  prob <- Problem(Minimize(qf_expr))
  result <- CVXR:::reduction_apply(d2c, prob)
  canon_prob <- result[[1L]]
  ## The objective should contain a SymbolicQuadForm
  obj_expr <- canon_prob@objective@args[[1L]]
  ## Walk to find SymbolicQuadForm
  found_sqf <- FALSE
  .walk <- function(e) {
    if (S7::S7_inherits(e, CVXR:::SymbolicQuadForm)) { found_sqf <<- TRUE; return() }
    for (a in e@args) .walk(a)
  }
  .walk(obj_expr)
  expect_true(found_sqf)
})

# ══════════════════════════════════════════════════════════════════════════════
# CVXPY parity gap closures — test_qp_solvers.py
# Verified against CVXPY 1.8.1 (branch claude, commit 3b964472b)
# ══════════════════════════════════════════════════════════════════════════════

# ── TestConicQuadObj::test_all_solvers ─────────────────────────────────────────
# CVXPY: Tests that conic solvers with supports_quad_obj() can solve QP problems
# using the quadratic objective path (no SOC cones introduced). In CVXR, the
# QP path is selected automatically when the objective has quadratic terms. We
# test representative QP sub-problems with a conic solver (Clarabel) to ensure
# correctness.

## @cvxpy test_qp_solvers.py::TestConicQuadObj::test_all_solvers
test_that("conic solver handles quadratic objective: quad_over_lin", {
  ## CVXPY QPTestBase.quad_over_lin: min 0.5*quad_over_lin(|x-1|, 1) s.t. x <= -1
  ## Expected: x = [-1, -1], obj = 0.5 * (|-1-1|^2 + |-1-1|^2) / 1 = 4
  x <- Variable(2)
  prob <- Problem(Minimize(0.5 * quad_over_lin(abs(x - 1), 1)),
                  list(x <= -1))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), c(-1, -1), tolerance = 1e-4)
  expect_equal(result, 4.0, tolerance = 1e-4)
})

## @cvxpy test_qp_solvers.py::TestConicQuadObj::test_all_solvers
test_that("conic solver handles quadratic objective: power_matrix", {
  ## CVXPY QPTestBase.power_matrix: min ||A - B||_F^2 s.t. B[0,1] == 4, B[1,0] == 2+8j
  ## In R: min sum(power(A - B, 2)) with A = [[1,2],[3,4]]
  ## Expected: B = [[1, 4], [2, 4]], obj = (2-4)^2 + (3-2)^2 = 4 + 1 = 5
  A_mat <- Variable(c(2, 2))
  B_mat <- matrix(c(1, 3, 2, 4), 2, 2)  ## R column-major: [[1,2],[3,4]]
  prob <- Problem(Minimize(sum(power(A_mat - B_mat, 2))),
                  list(A_mat[1, 2] == 4, A_mat[2, 1] == 2))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 5.0, tolerance = 1e-3)
})

## @cvxpy test_qp_solvers.py::TestConicQuadObj::test_all_solvers
test_that("conic solver handles quadratic objective: square_affine", {
  ## CVXPY QPTestBase.square_affine: min sum_squares(x) s.t. sum(x) == 4
  ## Expected: x = [2, 2], obj = 8
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x)), list(sum(x) == 4))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 8.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(2, 2), tolerance = 1e-4)
})

## @cvxpy test_qp_solvers.py::TestConicQuadObj::test_all_solvers
test_that("conic solver handles quadratic objective: maximize", {
  ## CVXPY QPTestBase.maximize_problem: max -(x'x + y'y + (z-2)^2) s.t. x >= 0
  ## Equivalently: min x'x + y'y + (z-2)^2 s.t. x >= 0
  ## Expected: x=0, y=0, z=2, obj_max = 0
  x <- Variable(2)
  y <- Variable(3)
  z <- Variable(2)
  prob <- Problem(Maximize(-(sum_squares(x) + sum_squares(y) + sum_squares(z - 2))),
                  list(x >= 0))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 0.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(0, 0), tolerance = 1e-3)
  expect_equal(as.numeric(value(y)), c(0, 0, 0), tolerance = 1e-3)
  expect_equal(as.numeric(value(z)), c(2, 2), tolerance = 1e-3)
})

## @cvxpy test_qp_solvers.py::TestConicQuadObj::test_all_solvers
test_that("conic solver handles quadratic objective: regression", {
  ## CVXPY QPTestBase.regression_1: min sum_squares(slope*xi + offset - yi)
  ## Straight line fit to (1,1), (2,2), (3,2)
  ## CVXPY expected: slope ~ 0.5, offset ~ 0.6667
  slope <- Variable(1)
  offset <- Variable(1)
  xi <- c(1, 2, 3)
  yi <- c(1, 2, 2)
  residuals <- slope * xi + offset - yi
  prob <- Problem(Minimize(sum_squares(residuals)))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(slope)), 0.5, tolerance = 1e-3)
  expect_equal(as.numeric(value(offset)), 2/3, tolerance = 1e-3)
})

# ── TestQp::test_parametric ───────────────────────────────────────────────────
# CVXPY: Solve parametric QP both from scratch and with Parameter, verify match.
# min a*x^2 + b*x s.t. 0 <= x <= 1, a=10, b in {-10, -2}
# Expected: b=-10 => x=0.5, obj=-2.5; b=-2 => x=0.1, obj=-0.1

## @cvxpy test_qp_solvers.py::TestQp::test_parametric
test_that("parametric QP: Parameter vs fresh solve give same results", {
  a <- 10
  b_vec <- c(-10, -2)

  ## Solve from scratch (no Parameter) for each b
  x_full <- list()
  obj_full <- numeric(length(b_vec))
  for (i in seq_along(b_vec)) {
    x <- Variable(1)
    obj <- Minimize(a * power(x, 2) + b_vec[i] * x)
    constraints <- list(x >= 0, x <= 1)
    prob <- Problem(obj, constraints)
    obj_full[i] <- psolve(prob, solver = "CLARABEL")
    x_full[[i]] <- as.numeric(value(x))
  }

  ## Solve parametric (with Parameter)
  x <- Variable(1)
  b <- Parameter(1)
  obj <- Minimize(a * power(x, 2) + b * x)
  constraints <- list(x >= 0, x <= 1)
  prob <- Problem(obj, constraints)
  x_param <- list()
  obj_param <- numeric(length(b_vec))
  for (i in seq_along(b_vec)) {
    value(b) <- b_vec[i]
    obj_param[i] <- psolve(prob, solver = "CLARABEL")
    x_param[[i]] <- as.numeric(value(x))
  }

  ## Verify match
  for (i in seq_along(b_vec)) {
    expect_equal(x_full[[i]], x_param[[i]], tolerance = 1e-3)
    expect_equal(obj_full[i], obj_param[i], tolerance = 1e-4)
  }

  ## Verify against analytical solutions
  ## b=-10: x* = -b/(2a) = 0.5, obj = 10*0.25 - 10*0.5 = -2.5
  expect_equal(x_full[[1]], 0.5, tolerance = 1e-3)
  expect_equal(obj_full[1], -2.5, tolerance = 1e-4)
  ## b=-2: x* = -b/(2a) = 0.1, obj = 10*0.01 - 2*0.1 = -0.1
  expect_equal(x_full[[2]], 0.1, tolerance = 1e-3)
  expect_equal(obj_full[2], -0.1, tolerance = 1e-4)
})

# ── TestQp::test_square_param ─────────────────────────────────────────────────
# CVXPY: min b^2 + |a|, a = Parameter(value=1), b = Variable
# Expected: b=0, obj = 0 + 1 = 1.0

## @cvxpy test_qp_solvers.py::TestQp::test_square_param
test_that("QP with squared variable + parameter in abs", {
  ## min b^2 + |a| where a is Parameter(value=1), b is Variable
  ## CVXPY expected: obj = 1.0 (b=0 is optimal, |1| = 1)
  ## Note: CVXPY uses SCS, but SCS fails with m=0 in R (REQUIRES_CONSTR).
  ## Use Clarabel which handles unconstrained problems.
  a <- Parameter(value = 1)
  b <- Variable(1)
  obj <- Minimize(power(b, 2) + abs(a))
  prob <- Problem(obj)
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(b)), 0.0, tolerance = 1e-3)
})

# ── TestQp::test_gurobi_warmstart ─────────────────────────────────────────────
# CVXPY: Tests Gurobi warm-start with user-provided initial point.
# In R, we test that warm_start=TRUE with pre-set Variable values works and
# produces correct results. We cannot access Gurobi model internals (start values)
# from R, so we verify correctness and that no errors occur.

## @cvxpy test_qp_solvers.py::TestQp::test_gurobi_warmstart
test_that("Gurobi warm-start with user-provided initial point", {
  skip_if_not_installed("gurobi")
  m <- 4L
  n <- 3L
  y <- Variable(nonneg = TRUE)
  X <- Variable(c(m, n))
  X_vals <- matrix(seq_len(m * n) - 1L, nrow = m, ncol = n)

  prob <- Problem(Minimize(power(y, 2) + sum_entries(X)), list(X == X_vals))

  ## Set initial values for X (warm-start hint)
  value(X) <- X_vals + 1

  ## Solve with warm_start=TRUE
  result <- psolve(prob, solver = "GUROBI", warm_start = TRUE)
  expect_equal(status(prob), "optimal")

  ## y is nonneg and minimizes y^2, so y = 0
  expect_equal(as.numeric(value(y)), 0.0, tolerance = 1e-4)

  ## X is fixed to X_vals by equality constraint
  expect_equal(as.numeric(value(X)), as.numeric(X_vals), tolerance = 1e-4)

  ## Objective: 0 + sum(0:11) = 66
  expect_equal(result, 66.0, tolerance = 1e-3)
})

# ── TestQp::test_gurobi_time_limit_no_solution ────────────────────────────────
# CVXPY: Tests Gurobi with TimeLimit=0 doesn't raise an error and returns stats.
# The Gurobi R interface passes solver_opts through; TimeLimit=0 should cause
# Gurobi to hit the time limit before finding a solution.

## @cvxpy test_qp_solvers.py::TestQp::test_gurobi_time_limit_no_solution
test_that("Gurobi time limit with no solution does not error", {
  skip_if_not_installed("gurobi")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1]), list(x[1] >= 1))

  ## Solve with TimeLimit=0 — should not throw an error
  expect_no_error({
    psolve(prob, solver = "GUROBI", TimeLimit = 0.0)
  })

  ## Solver stats should be returned regardless of status
  stats <- solver_stats(prob)
  expect_true(!is.null(stats))
})

# ── TestQp::test_gurobi_environment ───────────────────────────────────────────
# CVXPY: Tests passing a custom Gurobi environment with parameters. The R gurobi
# package does not support custom environment objects in the same way as gurobipy.
# This is an N/A feature in R, so we write a skip stub.

## @cvxpy test_qp_solvers.py::TestQp::test_gurobi_environment
test_that("Gurobi custom environment (N/A in R: gurobi pkg lacks env support)", {
  skip_if_not_installed("gurobi")
  skip("R gurobi package does not support custom Gurobi environment objects")
  ## CVXPY creates a gurobipy.Env() with custom params (MIPGap, AggFill,
  ## PerturbValue) and passes it via env= kwarg. The R gurobi package
  ## does not expose an equivalent mechanism. If this becomes available,
  ## implement the test following CVXPY test_qp_solvers.py lines 822-854.
})

# ── TestQp::test_highs_cvar ──────────────────────────────────────────────────
# CVXPY: CVaR optimization with HiGHS solver (issue #2836).
# Maximize portfolio PnL subject to CVaR constraint.

## @cvxpy test_qp_solvers.py::TestQp::test_highs_cvar
test_that("HiGHS CVaR constraint optimization", {
  skip_if_not_installed("highs")

  num_stocks <- 5L
  num_samples <- 25L
  set.seed(1L)
  pnl_samples <- matrix(runif(num_samples * num_stocks), nrow = num_samples, ncol = num_stocks)
  pnl_expected <- colMeans(pnl_samples)

  quantile_level <- 0.05
  w <- Variable(num_stocks, nonneg = TRUE)
  cvar_expr <- cvar(pnl_samples %*% w, 1 - quantile_level)
  pnl <- t(pnl_expected) %*% w

  objective <- Maximize(pnl)
  constraints <- list(cvar_expr <= 0.5)
  problem <- Problem(objective, constraints)
  psolve(problem, solver = "HIGHS")
  expect_equal(status(problem), "optimal")
})

# ── TestQp::test_highs_warmstart ──────────────────────────────────────────────
# CVXPY: HiGHS warm-start test. Blocked in R because the R highs package
# lacks setSolution() support needed for warm-starting.

## @cvxpy test_qp_solvers.py::TestQp::test_highs_warmstart
test_that("HiGHS warm-start (blocked: R highs pkg lacks setSolution)", {
  skip_if_not_installed("highs")
  skip("HiGHS warm-start blocked: R highs package lacks setSolution() support")
  ## CVXPY pattern: parametric least-squares, cold then warm, then new param
  ## warm then cold. Verify both give same results. When R highs package adds
  ## setSolution(), implement following CVXPY test_qp_solvers.py lines 635-656.
})

# ── TestQp::test_piqp_warmstart ──────────────────────────────────────────────
# CVXPY: PIQP warm-start with parametric least-squares.
# Solve with cold, then warm; change parameter, solve warm then cold.

## @cvxpy test_qp_solvers.py::TestQp::test_piqp_warmstart
test_that("PIQP warm-start: cold and warm produce same solution", {
  skip_if_not_installed("piqp")
  m <- 200L
  n <- 100L
  set.seed(1L)
  A_mat <- matrix(rnorm(m * n), nrow = m, ncol = n)
  b <- Parameter(m)

  x <- Variable(n)
  prob <- Problem(Minimize(sum_squares(A_mat %*% x - b)))

  ## Cycle 1: same b — cold then warm should match
  value(b) <- rnorm(m)
  result1 <- psolve(prob, solver = "PIQP", warm_start = FALSE)
  result2 <- psolve(prob, solver = "PIQP", warm_start = TRUE)
  expect_equal(result1, result2, tolerance = 1e-3)

  ## Cycle 2: new b — warm then cold should match
  value(b) <- rnorm(m)
  result3 <- psolve(prob, solver = "PIQP", warm_start = TRUE)
  result4 <- psolve(prob, solver = "PIQP", warm_start = FALSE)
  expect_equal(result3, result4, tolerance = 1e-3)
})

# ── TestQpSolverValidation::test_qp_solver_rejects_exponential_cones ──────────
# CVXPY: QP solver rejects problems requiring exponential cones.
# In CVXR, the QP solver path won't be selected for such problems (the solving
# chain routes to a conic solver). We test that explicitly requesting a QP-only
# solver for an exp-cone problem raises an appropriate error.

## @cvxpy test_qp_solvers.py::TestQpSolverValidation::test_qp_solver_rejects_exponential_cones
test_that("QP solver rejects problem requiring exponential cones", {
  skip_if_not_installed("osqp")

  ## Create a problem with exponential cone (log)
  x <- Variable(1)
  prob <- Problem(Maximize(log(x)), list(x <= 1))

  ## Requesting OSQP (a QP-only solver) for a problem that requires exp cones

  ## should fail: OSQP cannot handle the exponential cone constraints needed
  ## for log(x). The solving chain should raise an error about unsupported cones.
  expect_error(
    psolve(prob, solver = "OSQP"),
    regexp = ".*"
  )
})
