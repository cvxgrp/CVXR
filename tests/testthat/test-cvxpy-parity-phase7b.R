## CVXPY Parity Tests — Phase 7b
## =============================
## Ports tests from CVXPY's test_constraints.py, test_conic_solvers.py
## (solver_test_helpers.py), and test_problem.py.
##
## This file covers constraint semantics, solver benchmarks, and problem
## queries not already in test-cvxpy-solve-parity.R.
##
## CVXPY source: /Users/naras/GitHub/cvxpy branch claude

# ═══════════════════════════════════════════════════════════════════════
# PART 1: Constraint value / violation / residual
# Source: test_constraints.py
# ═══════════════════════════════════════════════════════════════════════

# ── Equality ──────────────────────────────────────────────────────────
## CVXPY SOURCE: test_constraints.py::test_equality()

## @cvxpy test_constraints.py::TestConstraints::test_equality
test_that("equality: shape and name", {
  x <- Variable(2)
  z <- Variable(2)
  constr <- x == z
  expect_equal(constr@shape, c(2L, 1L))
})

## @cvxpy test_constraints.py::TestConstraints::test_equality
test_that("equality: value TRUE when satisfied", {
  x <- Variable(2)
  z <- Variable(2)
  constr <- x == z
  save_leaf_value(x, matrix(c(2, 1), 2, 1))
  save_leaf_value(z, matrix(c(2, 1), 2, 1))
  expect_true(value(constr))
  expect_equal(as.numeric(violation(constr)), c(0, 0), tolerance = 1e-8)
})

## @cvxpy test_constraints.py::TestConstraints::test_equality
test_that("equality: value FALSE and violation when not satisfied", {
  x <- Variable(2)
  z <- Variable(2)
  constr <- x == z
  save_leaf_value(x, matrix(c(2, 1), 2, 1))
  save_leaf_value(z, matrix(c(2, 2), 2, 1))
  expect_false(value(constr))
  expect_equal(as.numeric(violation(constr)), c(0, 1), tolerance = 1e-8)
})

## @cvxpy test_constraints.py::TestConstraints::test_equality
test_that("equality: incompatible dimensions error", {
  x <- Variable(2)
  y <- Variable(3)
  expect_error(x == y)
})

# ── Inequality ────────────────────────────────────────────────────────
## CVXPY SOURCE: test_constraints.py::test_inequality()

## @cvxpy test_constraints.py::TestConstraints::test_inequality
test_that("inequality: shape and name", {
  x <- Variable(2)
  z <- Variable(2)
  constr <- x <= z
  expect_equal(constr@shape, c(2L, 1L))
})

## @cvxpy test_constraints.py::TestConstraints::test_inequality
test_that("inequality: value TRUE when satisfied", {
  x <- Variable(2)
  z <- Variable(2)
  constr <- x <= z
  save_leaf_value(x, matrix(c(1, 1), 2, 1))
  save_leaf_value(z, matrix(c(2, 2), 2, 1))
  expect_true(value(constr))
  expect_equal(as.numeric(violation(constr)), c(0, 0), tolerance = 1e-8)
})

## @cvxpy test_constraints.py::TestConstraints::test_inequality
test_that("inequality: value FALSE and violation when violated", {
  x <- Variable(2)
  z <- Variable(2)
  constr <- x <= z
  save_leaf_value(x, matrix(c(2, 1), 2, 1))
  save_leaf_value(z, matrix(c(2, 0), 2, 1))
  expect_false(value(constr))
  expect_equal(as.numeric(violation(constr)), c(0, 1), tolerance = 1e-8)
})

## @cvxpy test_constraints.py::TestConstraints::test_inequality
test_that("inequality: incompatible dimensions error", {
  x <- Variable(2)
  y <- Variable(3)
  expect_error(x <= y)
})

# ── >= operator ───────────────────────────────────────────────────────
## CVXPY SOURCE: test_constraints.py::test_geq()

## @cvxpy test_constraints.py::TestConstraints::test_geq
test_that("geq: z >= x creates x <= z", {
  x <- Variable(2)
  z <- Variable(2)
  constr <- z >= x
  expect_equal(constr@shape, c(2L, 1L))
})

# ── PSD constraint ────────────────────────────────────────────────────
## CVXPY SOURCE: test_constraints.py::test_psd_constraint()

## @cvxpy test_constraints.py::TestConstraints::test_psd_constraint
test_that("PSD: value TRUE for PSD matrix", {
  A <- Variable(c(2, 2))
  B <- Variable(c(2, 2))
  ## A - B >> 0 → PSD(A - B)
  constr <- PSD(A - B)
  expect_equal(constr@shape, c(2L, 2L))
  ## Set A = [[2,-1],[1,2]], B = [[1,0],[0,1]]
  ## A - B = [[1,-1],[1,1]], eigenvalues = 1±i (not PSD in real sense)
  ## Better: A = [[3,0],[0,3]], B = [[1,0],[0,1]] → A-B = [[2,0],[0,2]] PSD
  save_leaf_value(A, matrix(c(3, 0, 0, 3), 2, 2))
  save_leaf_value(B, matrix(c(1, 0, 0, 1), 2, 2))
  expect_true(value(constr))
  expect_equal(as.numeric(violation(constr)), 0, tolerance = 1e-8)
})

## @cvxpy test_constraints.py::TestConstraints::test_psd_constraint
test_that("PSD: value FALSE for non-PSD matrix", {
  A <- Variable(c(2, 2))
  B <- Variable(c(2, 2))
  constr <- PSD(A - B)
  ## A - B = [[-1,0],[0,-1]], eigenvalues = -1, -1 → not PSD
  save_leaf_value(A, matrix(c(0, 0, 0, 0), 2, 2))
  save_leaf_value(B, matrix(c(1, 0, 0, 1), 2, 2))
  expect_false(value(constr))
  expect_equal(as.numeric(violation(constr)), 1, tolerance = 1e-8)
})

## @cvxpy test_constraints.py::TestConstraints::test_psd_constraint
test_that("PSD: non-square matrix error", {
  x <- Variable(2)
  expect_error(PSD(x), "square")
})

# ── SOC constraint ────────────────────────────────────────────────────
## CVXPY SOURCE: test_constraints.py::test_soc_constraint()

## @cvxpy test_constraints.py::TestConstraints::test_soc_constraint
test_that("SOC: basic construction", {
  x <- Variable(2)
  z <- Variable(2)
  a <- Variable(1)
  b <- Variable(1)
  exp_val <- x + z
  scalar_exp <- a + b
  constr <- SOC(scalar_exp, exp_val)
  ## Should have 1 cone of size 3 (1 for t + 2 for x)
  expect_equal(cone_sizes(constr), 3L)
})

## @cvxpy test_constraints.py::TestConstraints::test_soc_constraint_scalar
test_that("SOC: scalar X solve (CVXPY issue #3054)", {
  ## min c s.t. |x| <= c, x == 3
  c_var <- Variable(1)
  x <- Variable(1)
  p <- Problem(Minimize(c_var), list(SOC(c_var, x), x == 3))
  psolve(p, verbose = FALSE)
  expect_equal(as.numeric(value(c_var)), 3, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 3, tolerance = 1e-3)
})

# ── PowCone3D constraint ─────────────────────────────────────────────
## CVXPY SOURCE: test_constraints.py::test_pow3d_constraint()

## @cvxpy test_constraints.py::TestConstraints::test_pow3d_constraint
test_that("PowCone3D: feasible residual near zero", {
  set.seed(0)
  n <- 3L
  alpha <- 0.275
  x <- Variable(n)
  y <- Variable(n)
  z <- Variable(n)
  con <- PowCone3D(x, y, z, alpha)
  ## Set feasible: z0 = x0^alpha * y0^(1-alpha)
  x0 <- 0.1 + runif(n)
  y0 <- 0.1 + runif(n)
  z0 <- x0^alpha * y0^(1 - alpha)
  z0[2] <- -z0[2]  # one negative z is OK
  save_leaf_value(x, matrix(x0, n, 1))
  save_leaf_value(y, matrix(y0, n, 1))
  save_leaf_value(z, matrix(z0, n, 1))
  ## residual returns a vector (one per cone element); max should be near 0
  expect_true(max(residual(con)) <= 1e-6)
})

## @cvxpy test_constraints.py::TestConstraints::test_pow3d_constraint
test_that("PowCone3D: invalid alpha", {
  n <- 3L
  x <- Variable(n)
  y <- Variable(n)
  z <- Variable(n)
  expect_error(PowCone3D(x, y, z, 1.001))
  expect_error(PowCone3D(x, y, z, -0.00001))
})

# ── PowConeND constraint ─────────────────────────────────────────────
## CVXPY SOURCE: test_constraints.py::test_pownd_constraint()

## @cvxpy test_constraints.py::TestConstraints::test_pownd_constraint
test_that("PowConeND: entries must sum to one", {
  n <- 4L
  W <- Variable(n)
  z <- Variable(1)
  set.seed(0)
  alpha <- 0.5 + runif(n)
  alpha <- alpha / sum(alpha)
  ## Entries don't sum to one
  expect_error(PowConeND(W, z, alpha + 0.01))
})

# ── NonNeg constraint ─────────────────────────────────────────────────
## CVXPY SOURCE: test_constraints.py::test_nonneg()

## @cvxpy test_constraints.py::TestConstraints::test_nonneg
test_that("NonNeg constraint solve", {
  x <- Variable(3)
  cc <- 0:2
  p <- Problem(Minimize(sum(x)), list(NonNeg(x - cc)))
  psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(as.numeric(value(x)), cc, tolerance = 1e-4)
})

# ── NonPos constraint ─────────────────────────────────────────────────
## CVXPY SOURCE: test_constraints.py::test_nonpos()

## @cvxpy test_constraints.py::TestConstraints::test_nonpos
test_that("NonPos constraint solve", {
  x <- Variable(3)
  cc <- 0:2
  p <- Problem(Maximize(sum(x)), list(NonPos(x - cc)))
  psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(as.numeric(value(x)), cc, tolerance = 1e-4)
})

# ── NonNeg dual ───────────────────────────────────────────────────────
## CVXPY SOURCE: test_constraints.py::test_nonneg_dual()

## @cvxpy test_constraints.py::TestConstraints::test_nonneg_dual
test_that("NonNeg dual matches Inequality dual", {
  x <- Variable(3)
  cc <- 0:2
  objective <- Minimize(sum(x))
  ## Reference with Inequality
  p1 <- Problem(objective, list(cc - x <= 0))
  psolve(p1, solver = CLARABEL_SOLVER, verbose = FALSE)
  dual_ref <- as.numeric(dual_value(p1@constraints[[1]]))
  ## Same with NonNeg
  p2 <- Problem(objective, list(NonNeg(x - cc)))
  psolve(p2, solver = CLARABEL_SOLVER, verbose = FALSE)
  dual_nn <- as.numeric(dual_value(p2@constraints[[1]]))
  expect_equal(dual_nn, dual_ref, tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════════
# PART 2: Standard solver benchmark problems
# Source: test_conic_solvers.py + solver_test_helpers.py
# ═══════════════════════════════════════════════════════════════════════

# ── lp_1: CVXOPT reference LP ─────────────────────────────────────────
## CVXPY SOURCE: solver_test_helpers.py::lp_1()
## min -4x[0] -5x[1] s.t. 2x[0]+x[1]<=3, x[0]+2x[1]<=3, x>=0
## opt = -9, x* = [1, 1], duals = [1, 2, 0, 0]

## @cvxpy NONE
test_that("Clarabel: lp_1 CVXOPT reference", {
  x <- Variable(2)
  objective <- Minimize(-4 * x[1] - 5 * x[2])
  constraints <- list(
    2 * x[1] + x[2] <= 3,
    x[1] + 2 * x[2] <= 3,
    x[1] >= 0,
    x[2] >= 0
  )
  p <- Problem(objective, constraints)
  result <- psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(result, -9, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-4)
  ## Dual values
  expect_equal(as.numeric(dual_value(constraints[[1]])), 1, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(constraints[[2]])), 2, tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(constraints[[3]])), 0, tolerance = 1e-2)
  expect_equal(as.numeric(dual_value(constraints[[4]])), 0, tolerance = 1e-2)
})

# ── lp_2: Bounded LP ─────────────────────────────────────────────────
## CVXPY SOURCE: solver_test_helpers.py::lp_2()
## min x[0]+0.5*x[1] s.t. x[0]>=-100, x[0]<=-10, x[1]==1
## opt = -99.5, x* = [-100, 1], duals = [1, 0, -0.5]

## @cvxpy NONE
test_that("Clarabel: lp_2 bounded LP", {
  x <- Variable(2)
  objective <- Minimize(x[1] + 0.5 * x[2])
  constraints <- list(
    x[1] >= -100,
    x[1] <= -10,
    x[2] == 1
  )
  p <- Problem(objective, constraints)
  result <- psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(result, -99.5, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(-100, 1), tolerance = 1e-3)
  expect_equal(as.numeric(dual_value(constraints[[1]])), 1, tolerance = 1e-2)
  expect_equal(as.numeric(dual_value(constraints[[3]])), -0.5, tolerance = 1e-2)
})

# ── lp_3: Unbounded LP ───────────────────────────────────────────────
## CVXPY SOURCE: solver_test_helpers.py::lp_3()
## min sum(x) s.t. x <= 1, x in R^5 → unbounded

## @cvxpy NONE
test_that("Clarabel: lp_3 unbounded", {
  x <- Variable(5)
  p <- Problem(Minimize(sum(x)), list(x <= 1))
  psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_true(status(p) %in% c(UNBOUNDED, INFEASIBLE_OR_UNBOUNDED))
})

## @cvxpy NONE
test_that("SCS: lp_3 unbounded", {
  x <- Variable(5)
  p <- Problem(Minimize(sum(x)), list(x <= 1))
  psolve(p, solver = SCS_SOLVER, verbose = FALSE)
  expect_true(status(p) %in% c(UNBOUNDED, INFEASIBLE_OR_UNBOUNDED))
})

# ── lp_4: Infeasible LP ──────────────────────────────────────────────
## CVXPY SOURCE: solver_test_helpers.py::lp_4()
## min sum(x) s.t. x <= 0, x >= 1 → infeasible

## @cvxpy NONE
test_that("Clarabel: lp_4 infeasible", {
  x <- Variable(5)
  p <- Problem(Minimize(sum(x)), list(x <= 0, x >= 1))
  psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(status(p), INFEASIBLE)
})

## @cvxpy NONE
test_that("SCS: lp_4 infeasible", {
  x <- Variable(5)
  p <- Problem(Minimize(sum(x)), list(x <= 0, x >= 1))
  psolve(p, solver = SCS_SOLVER, verbose = FALSE)
  expect_equal(status(p), INFEASIBLE)
})

# ── socp_0: Simple SOCP ──────────────────────────────────────────────
## CVXPY SOURCE: solver_test_helpers.py::socp_0()
## min ||x||_2 + 1 s.t. x == 0 → opt = 1, x* = [0, 0]

## @cvxpy NONE
test_that("Clarabel: socp_0", {
  x <- Variable(2)
  p <- Problem(Minimize(p_norm(x, 2) + 1), list(x == 0))
  result <- psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(result, 1, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(0, 0), tolerance = 1e-4)
})

# ── socp_1: SOCP with SOC constraint ─────────────────────────────────
## CVXPY SOURCE: solver_test_helpers.py::socp_1()
## min 3x[0]+2x[1]+x[2] s.t. ||x||_2<=y, x[0]+x[1]+3x[2]>=1, y<=5
## opt = -13.5486, x* = [-3.8746, -2.1298, 2.3348], y* = 5

## @cvxpy NONE
test_that("Clarabel: socp_1 with SOC constraint", {
  x <- Variable(3)
  y <- Variable(1)
  soc <- SOC(y, x)
  constraints <- list(soc,
                      x[1] + x[2] + 3 * x[3] >= 1.0,
                      y <= 5)
  obj <- Minimize(3 * x[1] + 2 * x[2] + x[3])
  p <- Problem(obj, constraints)
  result <- psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(result, -13.548638904065102, tolerance = 1e-3)
  xv <- as.numeric(value(x))
  expect_equal(xv[1], -3.87462, tolerance = 1e-3)
  expect_equal(xv[2], -2.12979, tolerance = 1e-3)
  expect_equal(xv[3], 2.33480, tolerance = 1e-3)
  expect_equal(as.numeric(value(y)), 5, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("SCS: socp_1 with SOC constraint", {
  x <- Variable(3)
  y <- Variable(1)
  soc <- SOC(y, x)
  constraints <- list(soc,
                      x[1] + x[2] + 3 * x[3] >= 1.0,
                      y <= 5)
  obj <- Minimize(3 * x[1] + 2 * x[2] + x[3])
  p <- Problem(obj, constraints)
  result <- psolve(p, solver = SCS_SOLVER, verbose = FALSE)
  expect_equal(result, -13.548638904065102, tolerance = 1e-2)
})

# ── sdp_1: Boyd & Vandenberghe Ex 8.3 ────────────────────────────────
## CVXPY SOURCE: solver_test_helpers.py::sdp_1('min')
## Correlation matrix with bounds on off-diagonal entries.
## min rho[0,3] → opt = -0.39

## @cvxpy NONE
test_that("Clarabel: sdp_1 minimize", {
  rho <- Variable(c(4, 4), symmetric = TRUE)
  constraints <- list(
    0.6 <= rho[1, 2], rho[1, 2] <= 0.9,
    0.8 <= rho[1, 3], rho[1, 3] <= 0.9,
    0.5 <= rho[2, 4], rho[2, 4] <= 0.7,
    -0.8 <= rho[3, 4], rho[3, 4] <= -0.4,
    rho[1, 1] == 1, rho[2, 2] == 1, rho[3, 3] == 1, rho[4, 4] == 1,
    PSD(rho)
  )
  obj <- Minimize(rho[1, 4])
  p <- Problem(obj, constraints)
  result <- psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(result, -0.39, tolerance = 0.02)
})

## @cvxpy NONE
test_that("SCS: sdp_1 minimize", {
  rho <- Variable(c(4, 4), symmetric = TRUE)
  constraints <- list(
    0.6 <= rho[1, 2], rho[1, 2] <= 0.9,
    0.8 <= rho[1, 3], rho[1, 3] <= 0.9,
    0.5 <= rho[2, 4], rho[2, 4] <= 0.7,
    -0.8 <= rho[3, 4], rho[3, 4] <= -0.4,
    rho[1, 1] == 1, rho[2, 2] == 1, rho[3, 3] == 1, rho[4, 4] == 1,
    PSD(rho)
  )
  obj <- Minimize(rho[1, 4])
  p <- Problem(obj, constraints)
  result <- psolve(p, solver = SCS_SOLVER, verbose = FALSE)
  expect_equal(result, -0.39, tolerance = 0.02)
})

## @cvxpy NONE
test_that("Clarabel: sdp_1 maximize", {
  rho <- Variable(c(4, 4), symmetric = TRUE)
  constraints <- list(
    0.6 <= rho[1, 2], rho[1, 2] <= 0.9,
    0.8 <= rho[1, 3], rho[1, 3] <= 0.9,
    0.5 <= rho[2, 4], rho[2, 4] <= 0.7,
    -0.8 <= rho[3, 4], rho[3, 4] <= -0.4,
    rho[1, 1] == 1, rho[2, 2] == 1, rho[3, 3] == 1, rho[4, 4] == 1,
    PSD(rho)
  )
  obj <- Maximize(rho[1, 4])
  p <- Problem(obj, constraints)
  result <- psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(result, 0.23, tolerance = 0.02)
})

# ── sdp_2: MOSEK SDO2 example ────────────────────────────────────────
## CVXPY SOURCE: solver_test_helpers.py::sdp_2()
## Two PSD variables: X1 (2x2), X2 (4x4) with trace objectives.
## opt = 52.40127214

## @cvxpy NONE
test_that("Clarabel: sdp_2 MOSEK example", {
  X1 <- Variable(c(2, 2), symmetric = TRUE)
  X2 <- Variable(c(4, 4), symmetric = TRUE)
  C1 <- matrix(c(1, 0, 0, 6), 2, 2)
  A1 <- matrix(c(1, 1, 1, 2), 2, 2)
  C2 <- matrix(c(1, -3, 0, 0, -3, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 4, 4)
  A2 <- matrix(c(0, 1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3), 4, 4)
  b_val <- 23
  k_val <- -3
  constraints <- list(
    matrix_trace(A1 %*% X1) + matrix_trace(A2 %*% X2) == b_val,
    X2[1, 2] <= k_val,
    PSD(X1),
    PSD(X2)
  )
  obj_expr <- Minimize(matrix_trace(C1 %*% X1) + matrix_trace(C2 %*% X2))
  p <- Problem(obj_expr, constraints)
  result <- psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(result, 52.40127214, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("SCS: sdp_2 MOSEK example", {
  X1 <- Variable(c(2, 2), symmetric = TRUE)
  X2 <- Variable(c(4, 4), symmetric = TRUE)
  C1 <- matrix(c(1, 0, 0, 6), 2, 2)
  A1 <- matrix(c(1, 1, 1, 2), 2, 2)
  C2 <- matrix(c(1, -3, 0, 0, -3, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 4, 4)
  A2 <- matrix(c(0, 1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3), 4, 4)
  b_val <- 23
  k_val <- -3
  constraints <- list(
    matrix_trace(A1 %*% X1) + matrix_trace(A2 %*% X2) == b_val,
    X2[1, 2] <= k_val,
    PSD(X1),
    PSD(X2)
  )
  obj_expr <- Minimize(matrix_trace(C1 %*% X1) + matrix_trace(C2 %*% X2))
  p <- Problem(obj_expr, constraints)
  result <- psolve(p, solver = SCS_SOLVER, verbose = FALSE)
  expect_equal(result, 52.40127214, tolerance = 0.1)
})

# ── expcone_1: Basic exponential cone ─────────────────────────────────
## CVXPY SOURCE: solver_test_helpers.py::expcone_1()
## min 3x[0]+2x[1]+x[2] s.t. 0.1<=sum(x)<=1, x>=0,
##   ExpCone(x[2],x[1],x[0])
## opt = 0.23534820622420757

## @cvxpy NONE
test_that("Clarabel: expcone_1", {
  x <- Variable(c(3, 1))
  cone_con <- ExpCone(x[3], x[2], x[1])
  constraints <- list(
    sum(x) <= 1.0,
    sum(x) >= 0.1,
    x >= 0,
    cone_con
  )
  obj <- Minimize(3 * x[1] + 2 * x[2] + x[3])
  p <- Problem(obj, constraints)
  result <- psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(result, 0.23534820622420757, tolerance = 1e-3)
  xv <- as.numeric(value(x))
  expect_equal(xv[1], 0.05462721, tolerance = 1e-3)
  expect_equal(xv[2], 0.02609378, tolerance = 1e-3)
  expect_equal(xv[3], 0.01927901, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("SCS: expcone_1", {
  x <- Variable(c(3, 1))
  cone_con <- ExpCone(x[3], x[2], x[1])
  constraints <- list(
    sum(x) <= 1.0,
    sum(x) >= 0.1,
    x >= 0,
    cone_con
  )
  obj <- Minimize(3 * x[1] + 2 * x[2] + x[3])
  p <- Problem(obj, constraints)
  result <- psolve(p, solver = SCS_SOLVER, verbose = FALSE)
  expect_equal(result, 0.23534820622420757, tolerance = 1e-2)
})

# ── pcp_2: PowCone3D for univariate power ────────────────────────────
## CVXPY SOURCE: solver_test_helpers.py::pcp_2()
## max (x^0.2)*(y^0.8) + z^0.4 - x via PowCone3D
## opt = -1.8073406786220672

## @cvxpy NONE
test_that("Clarabel: pcp_2 PowCone3D univariate power", {
  x <- Variable(3)
  hypos <- Variable(2)
  objective <- Minimize(-sum(hypos) + x[1])
  ## In CVXPY (1D): hstack([x[0], x[2]]) → (2,)
  ## In R (2D): vstack(x[1], x[3]) → (2, 1) to match hypos shape (2, 1)
  arg1 <- vstack(x[1], x[3])
  arg2 <- vstack(x[2], 1.0)
  pc_con <- PowCone3D(arg1, arg2, hypos, c(0.2, 0.4))
  constraints <- list(
    x[1] + x[2] + 0.5 * x[3] == 2,
    pc_con
  )
  p <- Problem(objective, constraints)
  result <- psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(result, -1.8073406786220672, tolerance = 1e-3)
  xv <- as.numeric(value(x))
  expect_equal(xv[1], 0.06393515, tolerance = 1e-2)
  expect_equal(xv[2], 0.78320961, tolerance = 1e-2)
  expect_equal(xv[3], 2.30571048, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("SCS: pcp_2 PowCone3D univariate power", {
  x <- Variable(3)
  hypos <- Variable(2)
  objective <- Minimize(-sum(hypos) + x[1])
  arg1 <- vstack(x[1], x[3])
  arg2 <- vstack(x[2], 1.0)
  pc_con <- PowCone3D(arg1, arg2, hypos, c(0.2, 0.4))
  constraints <- list(
    x[1] + x[2] + 0.5 * x[3] == 2,
    pc_con
  )
  p <- Problem(objective, constraints)
  result <- psolve(p, solver = SCS_SOLVER, verbose = FALSE)
  expect_equal(result, -1.8073406786220672, tolerance = 1e-2)
})

# ═══════════════════════════════════════════════════════════════════════
# PART 3: Problem queries
# Source: test_problem.py
# ═══════════════════════════════════════════════════════════════════════

# ── is_dcp ────────────────────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py::test_is_dcp()

## @cvxpy test_problem.py::TestProblem::test_is_dcp
test_that("is_dcp: minimize convex is DCP", {
  a <- Variable(1)
  p <- Problem(Minimize(norm_inf(a)))
  expect_true(is_dcp(p))
})

## @cvxpy test_problem.py::TestProblem::test_is_dcp
test_that("is_dcp: maximize convex is NOT DCP", {
  a <- Variable(1)
  p <- Problem(Maximize(norm_inf(a)))
  expect_false(is_dcp(p))
})

# ── is_lp ─────────────────────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py::test_is_lp()

## @cvxpy test_problem.py::TestProblem::test_is_lp
test_that("is_lp: linear objective + linear constraints", {
  set.seed(42)
  A <- matrix(rnorm(12), 4, 3)
  b <- rnorm(4)
  cc <- rnorm(3)
  Aeq <- matrix(rnorm(6), 2, 3)
  beq <- rnorm(2)
  y <- Variable(3)

  ## Simple LP
  p <- Problem(Minimize(t(cc) %*% y), list(A %*% y <= b))
  expect_true(is_lp(p))

  ## LP with equality
  p <- Problem(Minimize(t(cc) %*% y), list(A %*% y <= b, Aeq %*% y == beq))
  expect_true(is_lp(p))

  ## Maximize is also LP
  p <- Problem(Maximize(t(cc) %*% y), list(A %*% y <= b))
  expect_true(is_lp(p))

  ## QP is not LP (quadratic objective)
  p <- Problem(Minimize(sum_squares(y)), list(A %*% y <= b))
  expect_false(is_lp(p))

  ## SOC constraint makes it not LP
  t_var <- Variable(1)
  p <- Problem(Minimize(t(cc) %*% y), list(SOC(t_var, y)))
  expect_false(is_lp(p))

  ## ExpCone constraint makes it not LP
  p <- Problem(Minimize(t(cc) %*% y),
               list(ExpCone(y[1], y[2], y[3])))
  expect_false(is_lp(p))
})

## Fixed in Phase 9a (commit 7f811e9) — is_lp/is_qp now check PSD/NSD attributes.
## @cvxpy test_problem.py::TestProblem::test_is_lp
test_that("is_lp: PSD variable makes it not LP", {
  X <- Variable(c(2, 2), PSD = TRUE)
  p <- Problem(Minimize(matrix_trace(X)), list(X[1, 1] >= 1))
  expect_false(is_lp(p))
})

# ── is_qp ─────────────────────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py::test_is_qp()

## @cvxpy test_problem.py::TestProblem::test_is_qp
test_that("is_qp: quadratic objective + linear constraints", {
  set.seed(42)
  A <- matrix(rnorm(12), 4, 3)
  b <- rnorm(4)
  Aeq <- matrix(rnorm(6), 2, 3)
  beq <- rnorm(2)
  y <- Variable(3)

  ## sum_squares is QP
  obj <- sum_squares(A %*% y - b)
  p <- Problem(Minimize(obj), list())
  expect_true(is_qp(p))

  ## sum_squares + linear constraints is QP
  p <- Problem(Minimize(obj), list(Aeq %*% y == beq))
  expect_true(is_qp(p))
})

# ── variables ─────────────────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py::test_variables()

## @cvxpy test_problem.py::TestProblem::test_variables
test_that("variables: extracts all unique variables", {
  a <- Variable(1)
  x <- Variable(2)
  b <- Variable(1)
  A <- Variable(c(2, 2))
  p <- Problem(Minimize(a), list(a <= x, b <= A + 2))
  vars <- variables(p)
  ## Should contain a, x, b, A (4 variables)
  expect_equal(length(vars), 4L)
  ids <- vapply(vars, function(v) v@id, integer(1))
  expect_true(a@id %in% ids)
  expect_true(x@id %in% ids)
  expect_true(b@id %in% ids)
  expect_true(A@id %in% ids)
})

# ── constants ─────────────────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py::test_constants()

## @cvxpy test_problem.py::TestProblem::test_constants
test_that("constants: extracts all constants", {
  x <- Variable(2)
  c1 <- matrix(rnorm(2), 1, 2)
  c2 <- rnorm(2)
  p <- Problem(Minimize(c1 %*% x), list(x >= c2))
  consts <- constants(p)
  expect_true(length(consts) >= 2L)
})

# ── dual_variables (multi-constraint) ─────────────────────────────────
## CVXPY SOURCE: test_problem.py::test_dual_variables()

## @cvxpy test_problem.py::TestProblem::test_dual_variables
test_that("dual variables: LP with multiple constraints (Clarabel)", {
  x <- Variable(2)
  z <- Variable(2)
  ## R matrix(c(1,3,2,4), 2, 2) = [[1,2],[3,4]] (column-major)
  ## Matching CVXPY: np.array([[1,2],[3,4]]) @ z == [-1,-4]
  p <- Problem(
    Minimize(norm1(x + z)),
    list(
      x >= c(2, 3),
      matrix(c(1, 3, 2, 4), 2, 2) %*% z == c(-1, -4),
      p_norm(x + z, 2) <= 100
    )
  )
  result <- psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  ## Verified with CVXPY (uv run python): result=3.5, x=[2,3], z=[-2,0.5]
  expect_equal(result, 3.5, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), c(2, 3), tolerance = 1e-2)
  expect_equal(as.numeric(value(z)), c(-2, 0.5), tolerance = 1e-2)
  ## Dual values verified with CVXPY Clarabel
  expect_equal(as.numeric(dual_value(p@constraints[[1]]))[2], 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(dual_value(p@constraints[[3]])), 0, tolerance = 1e-2)
})

# ── indexing in solve ─────────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py::test_indexing()

## @cvxpy test_problem.py::TestProblem::test_indexing
test_that("indexing: vector variable element constraints", {
  x <- Variable(2)
  p <- Problem(Maximize(x[1]), list(x[1] <= 2, x[2] == 3))
  result <- psolve(p, verbose = FALSE)
  expect_equal(result, 2, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(2, 3), tolerance = 1e-4)
})

## @cvxpy test_problem.py::TestProblem::test_indexing
test_that("indexing: large matrix sum", {
  n <- 10L
  A_val <- matrix(seq_len(n * n) - 1L, n, n)  # 0..99
  X <- Variable(c(n, n))
  p <- Problem(Minimize(sum(X)), list(X == A_val))
  result <- psolve(p, verbose = FALSE)
  expected <- sum(0:(n * n - 1))
  expect_equal(result, expected, tolerance = 1e-2)
})

# ── scalar LP full test ──────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py::test_scalar_lp() (Test 3)

## @cvxpy test_problem.py::TestProblem::test_scalar_lp
test_that("scalar LP: three variables with offset (Clarabel)", {
  a <- Variable(1)
  b <- Variable(1)
  c_var <- Variable(1)
  p <- Problem(Minimize(3 * a - b + 100),
               list(a >= 2,
                    b + 5 * c_var - 2 == a,
                    b <= 5 + c_var))
  result <- psolve(p, solver = CLARABEL_SOLVER, verbose = FALSE)
  expect_equal(result, 101 + 1.0 / 6, tolerance = 1e-3)
  expect_equal(as.numeric(value(a)), 2, tolerance = 1e-3)
  expect_equal(as.numeric(value(b)), 5 - 1.0 / 6, tolerance = 1e-3)
  expect_equal(as.numeric(value(c_var)), -1.0 / 6, tolerance = 1e-3)
})

# ── Constant infeasible ──────────────────────────────────────────────
## CVXPY SOURCE: test_problem.py::test_constant_infeasible()

## @cvxpy test_problem.py::TestProblem::test_constant_infeasible
test_that("constant infeasible problem", {
  ## Zero-variable problem: ConstantSolver evaluates constraints directly
  p <- Problem(Maximize(0), list(Constant(0) == 1))
  psolve(p, verbose = FALSE)
  expect_equal(status(p), INFEASIBLE)
})

## @cvxpy test_problem.py::TestProblem::test_constant_infeasible
test_that("constant feasible problem", {
  ## Zero-variable problem where constraints are satisfied
  p <- Problem(Minimize(3), list(Constant(1) <= 2))
  result <- psolve(p, verbose = FALSE)
  expect_equal(status(p), OPTIMAL)
  expect_equal(result, 3)
})

# ── Solver stats ──────────────────────────────────────────────────────

## @cvxpy test_problem.py::TestProblem::test_solver_stats
test_that("solver_stats populated after solve", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(x)), list(x >= 1))
  psolve(p, verbose = FALSE, solver = "CLARABEL")
  stats <- solver_stats(p)
  expect_true(!is.null(stats))
  expect_true(!is.null(stats@solver_name))
  expect_true(stats@solve_time >= 0)
})

# ═══════════════════════════════════════════════════════════════════════
# PART 4: Power tools parity (validate R matches Python)
# Source: CVXPY cvxpy/utilities/power_tools.py
# ═══════════════════════════════════════════════════════════════════════

# ── limit_denominator (CPython algorithm) ─────────────────────────────

## @cvxpy test_power_tools.py::TestGeoMean::test_multi_step_dyad_completion
test_that("limit_denominator: matches CPython Fraction", {
  ## fractions.Fraction(7, 16).limit_denominator(5) == Fraction(2, 5)
  result <- .limit_denominator(as.bigq(7, 16), 5L)
  expect_equal(as.character(result), "2/5")

  ## fractions.Fraction(1, 3).limit_denominator(1) == Fraction(0, 1)
  result <- .limit_denominator(as.bigq(1, 3), 1L)
  expect_equal(as.character(result), "0")

  ## fractions.Fraction(1, 3).limit_denominator(2) == Fraction(1, 2)
  result <- .limit_denominator(as.bigq(1, 3), 2L)
  expect_equal(as.character(result), "1/2")

  ## Already within limit → unchanged
  result <- .limit_denominator(as.bigq(3, 7), 10L)
  expect_equal(as.character(result), "3/7")
})

# ── fracify ───────────────────────────────────────────────────────────
## CVXPY: fracify([1, 2, 3]) → (w, w_dyad)

## @cvxpy test_power_tools.py::TestGeoMean::test_multi_step_dyad_completion
test_that("fracify: integer weights [1,2,3]", {
  result <- fracify(c(1, 2, 3))
  w <- result$w
  w_dyad <- result$w_dyad
  ## w should sum to 1
  expect_equal(as.numeric(sum(w)), 1)
  ## w_dyad should sum to 1 and all denominators are powers of 2
  expect_equal(as.numeric(sum(w_dyad)), 1)
  expect_true(is_dyad_weight(w_dyad))
  ## Non-zero positions of w should be subset of non-zero positions of w_dyad
  expect_true(check_dyad(w, w_dyad))
})

## @cvxpy test_power_tools.py::TestGeoMean::test_multi_step_dyad_completion
test_that("fracify: equal weights [1,1,1,1,1]", {
  result <- fracify(c(1, 1, 1, 1, 1))
  w <- result$w
  w_dyad <- result$w_dyad
  expect_equal(as.numeric(sum(w)), 1)
  expect_true(is_dyad_weight(w_dyad))
  expect_true(check_dyad(w, w_dyad))
})

## @cvxpy test_power_tools.py::TestGeoMean::test_multi_step_dyad_completion
test_that("fracify: standard basis [0, 0, 1]", {
  result <- fracify(c(0, 0, 1))
  w <- result$w
  w_dyad <- result$w_dyad
  ## w = (0, 0, 1), w_dyad = (0, 0, 1) (already dyadic)
  expect_equal(as.numeric(w), c(0, 0, 1))
  expect_equal(as.numeric(w_dyad), c(0, 0, 1))
})

# ── dyad_completion ───────────────────────────────────────────────────

## @cvxpy test_power_tools.py::TestGeoMean::test_multi_step_dyad_completion
test_that("dyad_completion: (1/3, 1/3, 1/3) → 4 elements", {
  w <- as.bigq(c(1, 1, 1), c(3, 3, 3))
  w_dyad <- dyad_completion(w)
  expect_true(is_dyad_weight(w_dyad))
  expect_equal(length(w_dyad), 4L)
  ## (1/4, 1/4, 1/4, 1/4)
  expect_true(all(as.numeric(w_dyad) == 0.25))
})

## @cvxpy test_power_tools.py::TestGeoMean::test_multi_step_dyad_completion
test_that("dyad_completion: already dyadic → unchanged", {
  w <- as.bigq(c(1, 1), c(2, 2))
  w_dyad <- dyad_completion(w)
  expect_equal(length(w_dyad), 2L)
  expect_equal(as.character(w_dyad[1]), "1/2")
  expect_equal(as.character(w_dyad[2]), "1/2")
})

# ── decompose ─────────────────────────────────────────────────────────

## @cvxpy test_power_tools.py::TestGeoMean::test_multi_step_dyad_completion
test_that("decompose: simple dyadic weight", {
  w_dyad <- as.bigq(c(1, 1), c(2, 2))
  tree <- decompose(w_dyad)
  ## (1/2, 1/2) is a root that splits into two basis vectors
  expect_true(length(tree$keys) >= 1L)
})

# ── pow_high / pow_mid / pow_neg ──────────────────────────────────────

## @cvxpy NONE
test_that("pow_high: p=2, approx=TRUE", {
  result <- pow_high(2)
  ## CVXPY: 1/p = 1/2, so p_frac = 1/2, inv_p = 2 (integer)
  ## When inv_p is integer, returns int(inv_p) not bigq
  expect_equal(as.numeric(result$p), 2)
  w <- result$w
  expect_true(is.bigq(w))
  expect_equal(as.numeric(sum(w)), 1, tolerance = 1e-15)
})

## @cvxpy NONE
test_that("pow_high: p=2, approx=FALSE gives numeric", {
  result <- pow_high(2, approx = FALSE)
  expect_true(is.numeric(result$p) || is.integer(result$p))
  ## approx=FALSE returns original p
  expect_equal(as.numeric(result$p), 2)
  ## w = (1/p, 1 - 1/p) = (0.5, 0.5)
  expect_equal(as.numeric(result$w), c(0.5, 0.5))
})

## @cvxpy NONE
test_that("pow_mid: p=0.5, approx=TRUE gives bigq", {
  result <- pow_mid(0.5)
  expect_true(is.bigq(result$p))
  w <- result$w
  expect_true(is.bigq(w))
  expect_equal(as.numeric(sum(w)), 1, tolerance = 1e-15)
})

## @cvxpy NONE
test_that("pow_mid: p=0.5, approx=FALSE gives numeric", {
  result <- pow_mid(0.5, approx = FALSE)
  expect_true(is.numeric(result$p))
  expect_equal(result$p, 0.5)
  expect_equal(as.numeric(result$w), c(0.5, 0.5))
})

## @cvxpy NONE
test_that("pow_neg: p=-1, approx=TRUE gives bigq", {
  result <- pow_neg(-1)
  expect_true(is.bigq(result$p))
  w <- result$w
  expect_true(is.bigq(w))
  expect_equal(as.numeric(sum(w)), 1, tolerance = 1e-15)
})

## @cvxpy NONE
test_that("pow_neg: p=-1, approx=FALSE gives numeric", {
  result <- pow_neg(-1, approx = FALSE)
  expect_true(is.numeric(result$p))
  expect_equal(result$p, -1)
  ## w = (p/(p-1), -1/(p-1)) = (-1/-2, -1/-2) = (0.5, 0.5)
  expect_equal(as.numeric(result$w), c(0.5, 0.5))
})

# ── gm_constrs ────────────────────────────────────────────────────────

## @cvxpy test_power_tools.py::TestGeoMean::test_3d_power_cone_approx
test_that("gm_constrs: two equal weights produces 1 SOC constraint", {
  t_var <- Variable(c(1, 1))
  x1 <- Variable(c(1, 1))
  x2 <- Variable(c(1, 1))
  p_weights <- as.bigq(c(1, 1), c(2, 2))
  constrs <- gm_constrs(t_var, list(x1, x2), p_weights)
  ## (1/2, 1/2) → 1 SOC constraint (geometric mean of two variables)
  expect_equal(length(constrs), 1L)
  expect_true(S7_inherits(constrs[[1]], SOC))
})
