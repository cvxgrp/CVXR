## CVXPY Parity Tests: Kronecker product, quad_form, and Constraints
## ================================================================
## These tests mirror CVXPY test files to verify CVXR behavior matches:
##   - test_kron_canon.py (Kronecker product with Variable/Parameter/Constant)
##   - test_quad_form.py  (quadratic form edge cases)
##   - test_constraints.py (PowCone3D, PowConeND, chained constraints, bounds)
##
## Reference values obtained via `uv run python` against CVXPY 1.8.1.
## Random values are hardcoded (not relying on seed parity between R and Python).

library(testthat)
library(CVXR)

# =====================================================================
# KRON CANON TESTS (7 tests)
# =====================================================================
# CVXPY SOURCE: test_kron_canon.py
#
# CVXPY's make_kron_prob(z_dims, c_dims, param, var_left, seed):
#   np.random.seed(seed=0) -> generates C_value, then L.
#   Optimal Z = L (because kron(C, Z) >= kron(C, L) with min sum(Z)).
#
# Python np.random.seed(0) values (row-major):
#   C(1,1) = [[0.55]],              L(2,2) = [[0.72,0.60],[0.54,0.42]]
#   C(2,1) = [[0.55],[0.72]],       L(2,2) = [[0.60,0.54],[0.42,0.65]]
#   C(1,2) = [[0.55,0.72]],         L(2,2) = [[0.60,0.54],[0.42,0.65]]
#   C(2,2) = [[0.55,0.72],[0.60,0.54]], L(2,2) = [[0.42,0.65],[0.44,0.89]]

# ── 1. test_gen_kronr_param ──────────────────────────────────────────
# CVXPY: kron(Parameter C, Variable Z) >= kron(C, L), Minimize sum(Z)
# Result: Z = L for all c_dims

## @cvxpy test_kron_canon.py::TestKronRightVar::test_gen_kronr_param
test_that("kron: gen_kronr_param — kron(Parameter, Variable) across c_dims", {
  ## c_dims = (1,1): C=0.55, L=[[0.72,0.60],[0.54,0.42]]
  C_param <- Parameter(c(1, 1))
  value(C_param) <- matrix(0.55, 1, 1)
  L <- matrix(c(0.72, 0.54, 0.60, 0.42), 2, 2)  # col-major
  Z <- Variable(c(2, 2))
  prob <- Problem(Minimize(sum_entries(Z)),
                  list(kron(C_param, Z) >= kron(C_param, Constant(L)), Z >= 0))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, sum(L), tolerance = 1e-3)
  expect_equal(value(Z), L, tolerance = 1e-3)

  ## c_dims = (2,1): C=[[0.55],[0.72]], L=[[0.60,0.54],[0.42,0.65]]
  C_param2 <- Parameter(c(2, 1))
  value(C_param2) <- matrix(c(0.55, 0.72), 2, 1)
  L2 <- matrix(c(0.60, 0.42, 0.54, 0.65), 2, 2)
  Z2 <- Variable(c(2, 2))
  prob2 <- Problem(Minimize(sum_entries(Z2)),
                   list(kron(C_param2, Z2) >= kron(C_param2, Constant(L2)), Z2 >= 0))
  val2 <- psolve(prob2, solver = CLARABEL_SOLVER)
  expect_equal(status(prob2), OPTIMAL)
  expect_equal(val2, sum(L2), tolerance = 1e-3)
  expect_equal(value(Z2), L2, tolerance = 1e-3)

  ## c_dims = (1,2): C=[[0.55,0.72]], L=[[0.60,0.54],[0.42,0.65]]
  C_param3 <- Parameter(c(1, 2))
  value(C_param3) <- matrix(c(0.55, 0.72), 1, 2)
  L3 <- matrix(c(0.60, 0.42, 0.54, 0.65), 2, 2)
  Z3 <- Variable(c(2, 2))
  prob3 <- Problem(Minimize(sum_entries(Z3)),
                   list(kron(C_param3, Z3) >= kron(C_param3, Constant(L3)), Z3 >= 0))
  val3 <- psolve(prob3, solver = CLARABEL_SOLVER)
  expect_equal(status(prob3), OPTIMAL)
  expect_equal(val3, sum(L3), tolerance = 1e-3)
  expect_equal(value(Z3), L3, tolerance = 1e-3)

  ## c_dims = (2,2): C=[[0.55,0.72],[0.60,0.54]], L=[[0.42,0.65],[0.44,0.89]]
  C_param4 <- Parameter(c(2, 2))
  value(C_param4) <- matrix(c(0.55, 0.60, 0.72, 0.54), 2, 2)
  L4 <- matrix(c(0.42, 0.44, 0.65, 0.89), 2, 2)
  Z4 <- Variable(c(2, 2))
  prob4 <- Problem(Minimize(sum_entries(Z4)),
                   list(kron(C_param4, Z4) >= kron(C_param4, Constant(L4)), Z4 >= 0))
  val4 <- psolve(prob4, solver = CLARABEL_SOLVER)
  expect_equal(status(prob4), OPTIMAL)
  expect_equal(val4, sum(L4), tolerance = 1e-3)
  expect_equal(value(Z4), L4, tolerance = 1e-3)
})

# ── 2. test_gen_kronr_const ──────────────────────────────────────────
# CVXPY: kron(Constant C, Variable Z) >= kron(C, L), Minimize sum(Z)

## @cvxpy test_kron_canon.py::TestKronLeftVar::test_gen_kronr_const
test_that("kron: gen_kronr_const — kron(Constant, Variable) across c_dims", {
  ## c_dims = (1,1)
  C1 <- Constant(matrix(0.55, 1, 1))
  L1 <- matrix(c(0.72, 0.54, 0.60, 0.42), 2, 2)
  Z1 <- Variable(c(2, 2))
  prob1 <- Problem(Minimize(sum_entries(Z1)),
                   list(kron(C1, Z1) >= kron(C1, Constant(L1)), Z1 >= 0))
  val1 <- psolve(prob1, solver = CLARABEL_SOLVER)
  expect_equal(status(prob1), OPTIMAL)
  expect_equal(val1, sum(L1), tolerance = 1e-3)
  expect_equal(value(Z1), L1, tolerance = 1e-3)

  ## c_dims = (2,2)
  C4 <- Constant(matrix(c(0.55, 0.60, 0.72, 0.54), 2, 2))
  L4 <- matrix(c(0.42, 0.44, 0.65, 0.89), 2, 2)
  Z4 <- Variable(c(2, 2))
  prob4 <- Problem(Minimize(sum_entries(Z4)),
                   list(kron(C4, Z4) >= kron(C4, Constant(L4)), Z4 >= 0))
  val4 <- psolve(prob4, solver = CLARABEL_SOLVER)
  expect_equal(status(prob4), OPTIMAL)
  expect_equal(val4, sum(L4), tolerance = 1e-3)
  expect_equal(value(Z4), L4, tolerance = 1e-3)
})

# ── 3. test_symvar_kronl_param ───────────────────────────────────────
# CVXPY: kron(symmetric Variable X, Parameter b) with bounds L <= kron(X,b) <= U

## @cvxpy test_kron_canon.py::TestKronLeftVar::test_symvar_kronl_param
test_that("kron: symvar_kronl_param — kron(symmetric Var, Parameter)", {
  ## Minimize sum(X) with symmetric X
  X <- Variable(c(2, 2), symmetric = TRUE)
  b <- Parameter(c(1, 1))
  value(b) <- matrix(1.5, 1, 1)
  L <- matrix(c(0.5, 2, 1, 3), 2, 2)   # col-major: [[0.5,1],[2,3]]
  U <- matrix(c(10, 12, 11, 13), 2, 2)  # col-major: [[10,11],[12,13]]
  kronX <- kron(X, b)

  prob_min <- Problem(Minimize(sum_entries(X)), list(U >= kronX, kronX >= L))
  val_min <- psolve(prob_min, solver = CLARABEL_SOLVER)
  expect_equal(status(prob_min), OPTIMAL)
  ## X_min = [[0.5, 2], [2, 3]] / 1.5 = [[0.333, 1.333], [1.333, 2.0]]
  ## (symmetry forces off-diagonal = max(1, 2)/1.5 = 1.333)
  Xval_min <- value(X)
  expect_equal(Xval_min[1, 1], 0.5 / 1.5, tolerance = 1e-3)
  expect_equal(Xval_min[2, 2], 3.0 / 1.5, tolerance = 1e-3)
  expect_equal(Xval_min[1, 2], Xval_min[2, 1], tolerance = 1e-6)  # symmetric

  ## Maximize sum(X)
  X2 <- Variable(c(2, 2), symmetric = TRUE)
  b2 <- Parameter(c(1, 1))
  value(b2) <- matrix(1.5, 1, 1)
  kronX2 <- kron(X2, b2)
  prob_max <- Problem(Maximize(sum_entries(X2)), list(U >= kronX2, kronX2 >= L))
  val_max <- psolve(prob_max, solver = CLARABEL_SOLVER)
  expect_equal(status(prob_max), OPTIMAL)
  ## X_max = [[10, 11], [11, 13]] / 1.5 (symmetry forces off-diag to min(11,12)/1.5)
  Xval_max <- value(X2)
  expect_equal(Xval_max[1, 2], Xval_max[2, 1], tolerance = 1e-6)
  expect_equal(val_max, 30.0, tolerance = 1e-2)
})

# ── 4. test_symvar_kronl_const ───────────────────────────────────────
# CVXPY: kron(symmetric Variable X, Constant b)

## @cvxpy test_kron_canon.py::TestKronLeftVar::test_symvar_kronl_const
test_that("kron: symvar_kronl_const — kron(symmetric Var, Constant)", {
  X <- Variable(c(2, 2), symmetric = TRUE)
  b <- Constant(matrix(1.5, 1, 1))
  L <- matrix(c(0.5, 2, 1, 3), 2, 2)
  U <- matrix(c(10, 12, 11, 13), 2, 2)
  kronX <- kron(X, b)

  ## Minimize
  prob_min <- Problem(Minimize(sum_entries(X)), list(U >= kronX, kronX >= L))
  val_min <- psolve(prob_min, solver = CLARABEL_SOLVER)
  expect_equal(status(prob_min), OPTIMAL)
  expect_equal(val_min, 5.0, tolerance = 1e-2)
  Xval <- value(X)
  expect_equal(Xval[1, 2], Xval[2, 1], tolerance = 1e-6)

  ## Maximize
  X2 <- Variable(c(2, 2), symmetric = TRUE)
  kronX2 <- kron(X2, b)
  prob_max <- Problem(Maximize(sum_entries(X2)), list(U >= kronX2, kronX2 >= L))
  val_max <- psolve(prob_max, solver = CLARABEL_SOLVER)
  expect_equal(status(prob_max), OPTIMAL)
  expect_equal(val_max, 30.0, tolerance = 1e-2)
})

# ── 5. test_scalar_kronl_param ───────────────────────────────────────
# CVXPY: kron(scalar Variable y, Parameter A)

## @cvxpy test_kron_canon.py::TestKronLeftVar::test_scalar_kronl_param
test_that("kron: scalar_kronl_param — kron(scalar Var, Parameter)", {
  y <- Variable(c(1, 1))
  A <- Parameter(c(2, 2))
  A_val <- matrix(c(1, 3, 2, 4), 2, 2)  # col-major: [[1,2],[3,4]]
  value(A) <- A_val
  L <- matrix(c(0.5, 2, 1, 3), 2, 2)
  U <- matrix(c(10, 12, 11, 13), 2, 2)
  krony <- kron(y, A)

  ## Min y s.t. U >= kron(y, A) >= L
  prob_min <- Problem(Minimize(y), list(U >= krony, krony >= L))
  val_min <- psolve(prob_min, solver = CLARABEL_SOLVER)
  expect_equal(status(prob_min), OPTIMAL)
  ## y_min = max(L / A_val) = max(0.5/1, 2/3, 1/2, 3/4) = 0.75
  expect_equal(val_min, 0.75, tolerance = 1e-3)

  ## Max y
  y2 <- Variable(c(1, 1))
  krony2 <- kron(y2, A)
  prob_max <- Problem(Maximize(y2), list(U >= krony2, krony2 >= L))
  val_max <- psolve(prob_max, solver = CLARABEL_SOLVER)
  expect_equal(status(prob_max), OPTIMAL)
  ## y_max = min(U / A_val) = min(10/1, 12/3, 11/2, 13/4) = 3.25
  expect_equal(val_max, 3.25, tolerance = 1e-3)
})

# ── 6. test_scalar_kronl_const ───────────────────────────────────────
# CVXPY: kron(scalar Variable y, Constant A)

## @cvxpy test_kron_canon.py::TestKronLeftVar::test_scalar_kronl_const
test_that("kron: scalar_kronl_const — kron(scalar Var, Constant)", {
  y <- Variable(c(1, 1))
  A_val <- matrix(c(1, 3, 2, 4), 2, 2)
  A <- Constant(A_val)
  L <- matrix(c(0.5, 2, 1, 3), 2, 2)
  U <- matrix(c(10, 12, 11, 13), 2, 2)
  krony <- kron(y, A)

  prob_min <- Problem(Minimize(y), list(U >= krony, krony >= L))
  val_min <- psolve(prob_min, solver = CLARABEL_SOLVER)
  expect_equal(status(prob_min), OPTIMAL)
  expect_equal(val_min, 0.75, tolerance = 1e-3)

  y2 <- Variable(c(1, 1))
  krony2 <- kron(y2, A)
  prob_max <- Problem(Maximize(y2), list(U >= krony2, krony2 >= L))
  val_max <- psolve(prob_max, solver = CLARABEL_SOLVER)
  expect_equal(status(prob_max), OPTIMAL)
  expect_equal(val_max, 3.25, tolerance = 1e-3)
})

# ── 7. test_gen_kronl_param ──────────────────────────────────────────
# CVXPY: kron(Variable Z, Parameter C) >= kron(L, C), Minimize sum(Z)

## @cvxpy test_kron_canon.py::TestKronLeftVar::test_gen_kronl_param
test_that("kron: gen_kronl_param — kron(Variable, Parameter) across c_dims", {
  ## c_dims = (1,1): C=0.55, L=[[0.72,0.60],[0.54,0.42]]
  C_param <- Parameter(c(1, 1))
  value(C_param) <- matrix(0.55, 1, 1)
  L <- matrix(c(0.72, 0.54, 0.60, 0.42), 2, 2)
  Z <- Variable(c(2, 2))
  prob <- Problem(Minimize(sum_entries(Z)),
                  list(kron(Z, C_param) >= kron(Constant(L), C_param), Z >= 0))
  val <- psolve(prob, solver = CLARABEL_SOLVER)
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, sum(L), tolerance = 1e-3)
  expect_equal(value(Z), L, tolerance = 1e-3)

  ## c_dims = (2,2): C=[[0.55,0.72],[0.60,0.54]], L=[[0.42,0.65],[0.44,0.89]]
  C_param4 <- Parameter(c(2, 2))
  value(C_param4) <- matrix(c(0.55, 0.60, 0.72, 0.54), 2, 2)
  L4 <- matrix(c(0.42, 0.44, 0.65, 0.89), 2, 2)
  Z4 <- Variable(c(2, 2))
  prob4 <- Problem(Minimize(sum_entries(Z4)),
                   list(kron(Z4, C_param4) >= kron(Constant(L4), C_param4), Z4 >= 0))
  val4 <- psolve(prob4, solver = CLARABEL_SOLVER)
  expect_equal(status(prob4), OPTIMAL)
  expect_equal(val4, sum(L4), tolerance = 1e-3)
  expect_equal(value(Z4), L4, tolerance = 1e-3)
})


# =====================================================================
# QUAD FORM TESTS (5 tests)
# =====================================================================
# CVXPY SOURCE: test_quad_form.py

# ── 1. test_singular_quad_form ───────────────────────────────────────
# CVXPY: quad_form with rank-deficient P on simplex; optimal at null-space direction

## @cvxpy test_quad_form.py::TestNonOptimal::test_singular_quad_form
test_that("quad_form: singular (rank-deficient) P on simplex", {
  ## Uses Python np.random.seed(1234) values hardcoded:
  ## v = [0.26283272, 0.04985441, 0.68731287] (normalized exp(randn))
  ## Q = E @ (A @ A^T) @ E^T where E = I - v v^T / (v^T v), rank = n-1
  ## Optimal: x = v, value = 0

  set.seed(1234)
  n <- 3L
  v <- exp(rnorm(n))
  v <- v / sum(v)

  A <- matrix(rnorm(n * n), n, n)
  Q <- A %*% t(A)

  ## Project onto orth complement of v

  E <- diag(n) - outer(v, v) / sum(v * v)
  Q <- E %*% Q %*% t(E)
  ## Symmetrize (numerical reasons)
  Q <- (Q + t(Q)) / 2
  expect_equal(as.integer(Matrix::rankMatrix(Q)), n - 1L)

  x <- Variable(n)
  q <- quad_form(x, Q)
  prob <- Problem(Minimize(q), list(x >= 0, sum(x) == 1))
  psolve(prob, solver = "OSQP")
  xopt <- as.numeric(value(x))
  yopt <- as.numeric(t(xopt) %*% Q %*% xopt)
  expect_equal(yopt, 0, tolerance = 1e-3)
  expect_equal(xopt, v, tolerance = 1e-3)
})

# ── 2. test_param_quad_form ──────────────────────────────────────────
# CVXPY: quad_form(x, P) where P is a PSD Parameter

## @cvxpy test_quad_form.py::TestNonOptimal::test_param_quad_form
test_that("quad_form: Parameter as P matrix (PSD)", {
  P <- Parameter(c(2, 2), PSD = TRUE)
  x <- Variable(2)
  cost <- quad_form(x, P)
  value(P) <- diag(2)
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  ## SCS matches CVXPY (which uses SCS for param quad_form)
  suppressWarnings(val <- psolve(prob, solver = "SCS"))
  ## x^T I x = 1 + 4 = 5
  expect_equal(val, 5, tolerance = 1e-1)
})

# ── 3. test_non_symmetric ────────────────────────────────────────────
# CVXPY: non-symmetric P raises error

## @cvxpy test_quad_form.py::TestNonOptimal::test_non_symmetric
test_that("quad_form: non-symmetric P raises error", {
  P <- matrix(c(2, 2, 3, 4), 2, 2)  # col-major: [[2,3],[2,4]] (not symmetric)
  x <- Variable(2)
  expect_error(quad_form(x, P), "symmetric|Hermitian")
})

# ── 4. test_non_psd ─────────────────────────────────────────────────
# CVXPY: indefinite P: quad_form constructs but Minimize fails DCP

## @cvxpy test_quad_form.py::TestNonOptimal::test_non_psd
test_that("quad_form: indefinite P is not DCP for Minimize", {
  P <- matrix(c(1, 0, 0, -1), 2, 2)
  x <- Variable(2)
  ## Forming quad_form with indefinite P warns but succeeds
  suppressWarnings({
    cost <- quad_form(x, P)
  })
  ## Not convex (indefinite)
  expect_false(is_convex(cost))
  ## Not concave either
  expect_false(is_concave(cost))
  ## Minimize should fail DCP
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  expect_error(psolve(prob, solver = "SCS"), "DCP")
})

# ── 5. test_assume_psd ──────────────────────────────────────────────
# CVXPY: quad_form(x, P, assume_PSD=True) forces convexity even if P is not PSD
# NOTE: CVXR does not yet support assume_PSD. This test documents the gap.

## @cvxpy test_quad_form.py::TestNonOptimal::test_assume_psd
test_that("quad_form: assume_PSD flag behavior (gap documentation)", {
  x <- Variable(3)
  A <- diag(3)
  ## Standard PSD: always convex
  expr <- quad_form(x, A)
  expect_true(is_convex(expr))

  ## NSD matrix: quad_form is concave, not convex
  A_neg <- -diag(3)
  suppressWarnings({
    expr_neg <- quad_form(x, A_neg)
  })
  expect_true(is_concave(expr_neg))
  expect_false(is_convex(expr_neg))

  ## If assume_PSD existed in CVXR, quad_form(x, -I, assume_PSD=TRUE)
  ## would report is_convex()=TRUE. Document this for future implementation.
  ## See CVXPY test_quad_form.py::test_assume_psd
})


# =====================================================================
# CONSTRAINT TESTS (5 tests)
# =====================================================================
# CVXPY SOURCE: test_constraints.py

# ── 1. test_pow3d_constraint ─────────────────────────────────────────
# CVXPY: PowCone3D construction, residual computation, and alpha validation

## @cvxpy test_constraints.py::TestConstraints::test_pow3d_constraint
test_that("PowCone3D: construction, residual, and solve (CVXPY parity)", {
  n <- 3L
  alpha <- 0.275

  x <- Variable(c(n, 1L))
  y <- Variable(c(n, 1L))
  z <- Variable(c(n, 1L))
  con <- PowCone3D(x, y, z, alpha)

  ## Feasible values (from np.random.seed(0)):
  ## x0 = [0.6488, 0.8152, 0.7028]
  ## y0 = [0.6449, 0.5237, 0.7459]
  ## z0 = x0^alpha * y0^(1-alpha), then z0[2] *= -1
  x0 <- c(0.6488135, 0.81518937, 0.70276338)
  y0 <- c(0.64488318, 0.5236548, 0.74589411)
  z0 <- x0^alpha * y0^(1 - alpha)
  z0[2] <- -z0[2]

  value(x) <- matrix(x0, n, 1)
  value(y) <- matrix(y0, n, 1)
  value(z) <- matrix(z0, n, 1)
  viol <- residual(con)
  expect_true(all(viol <= 1e-7))

  ## Infeasible: make x[1] negative
  x1 <- x0
  x1[1] <- -0.9 * x1[1]
  value(x) <- matrix(x1, n, 1)
  viol2 <- residual(con)
  expect_true(max(viol2) >= 0.99 * abs(x1[1]))

  ## Invalid alpha values
  expect_error(PowCone3D(x, y, z, 1.001))
  expect_error(PowCone3D(x, y, z, -0.00001))

  ## Solve test: PowCone3D scalar alpha broadcast
  ## CVXPY: minimize ||x - x0|| s.t. PowCone3D(x0[1], x0[2], x0[3], 0.25), x <= -10
  ## Expected value: 17.3205 (approx 10*sqrt(3))
  x0_var <- Variable(3)
  x_var <- Variable(3)
  cons <- list(PowCone3D(x0_var[1, 1], x0_var[2, 1], x0_var[3, 1], 0.25),
               x_var <= -10)
  obj <- Minimize(p_norm(x_var - x0_var))
  prob <- Problem(obj, cons)
  result <- psolve(prob)
  expect_equal(result, 17.320508, tolerance = 1e-2)
})

# ── 2. test_pownd_constraint ─────────────────────────────────────────
# CVXPY: PowConeND validation and violation computation

## @cvxpy test_constraints.py::TestConstraints::test_pownd_constraint
test_that("PowConeND: validation errors and violation (CVXPY parity)", {
  n <- 4L
  W <- Variable(c(n, 1L))
  z <- Variable(1L)

  ## Alpha from np.random.seed(0), normalized:
  alpha <- c(0.23773727, 0.27545012, 0.24996623, 0.23684638)

  ## Error: entries don't sum to 1
  bad_alpha <- alpha + 0.01
  expect_error(
    PowConeND(W, z, Constant(matrix(bad_alpha, n, 1)), axis = 2L),
    "sum to 1"
  )

  ## Error: shapes don't match (row vs column)
  expect_error(
    PowConeND(W, z, Constant(matrix(alpha, 1, n)), axis = 2L)
  )

  ## Compute violation
  con <- PowConeND(W, z, Constant(matrix(alpha, n, 1)), axis = 2L)
  W0 <- c(0.5236548, 0.74589411, 0.53758721, 0.991773)
  z0 <- prod(W0^alpha) + 0.05  # slightly infeasible
  value(W) <- matrix(W0, n, 1)
  value(z) <- z0
  viol <- violation(con)
  ## CVXPY: violation in [0.01, 0.06]
  expect_true(viol >= 0.01)
  expect_true(viol <= 0.06)
})

# ── 3. test_chained_constraints ──────────────────────────────────────
# CVXPY: chaining constraints (z <= x <= 1) should error

## @cvxpy test_constraints.py::TestConstraints::test_chained_constraints
test_that("Chained constraints raise error (CVXPY parity)", {
  x <- Variable(2, name = "x")
  z <- Variable(2, name = "z")

  ## In R, z <= x produces a Constraint; applying <= 1 to it should error
  constr <- (z <= x)
  expect_error(constr <= 1)

  ## Similarly for equality chaining
  constr_eq <- (x == z)
  expect_error(constr_eq == 1)
})

# ── 4. test_bound_properties ─────────────────────────────────────────
# CVXPY: Variable bounds attribute (lower, upper, both)

## @cvxpy test_constraints.py::TestConstraints::test_bound_properties
test_that("Variable bounds properties (CVXPY parity)", {
  ## Lower bound only
  x1 <- Variable(c(1, 1), bounds = list(1, NULL))
  attrs1 <- x1@attributes
  expect_true(all(attrs1$bounds[[1L]] > -Inf))
  expect_true(all(attrs1$bounds[[2L]] == Inf))

  ## Upper bound only
  x2 <- Variable(c(1, 1), bounds = list(NULL, 1))
  attrs2 <- x2@attributes
  expect_true(all(attrs2$bounds[[1L]] == -Inf))
  expect_true(all(attrs2$bounds[[2L]] < Inf))

  ## Both bounds
  x3 <- Variable(c(1, 1), bounds = list(1, 2))
  attrs3 <- x3@attributes
  expect_equal(attrs3$bounds[[1L]], 1)
  expect_equal(attrs3$bounds[[2L]], 2)
})

# ── 5. test_bounds_attr ──────────────────────────────────────────────
# CVXPY: bounds as attribute vs explicit constraints give same results

## @cvxpy test_constraints.py::TestConstraints::test_bounds_attr
test_that("Variable bounds as attribute vs explicit constraints (CVXPY parity)", {
  ## Bounded variable: minimize sum(x) with bounds [0, 5]
  x_bounded <- Variable(2, bounds = list(0, 5))
  prob_bounded <- Problem(Minimize(sum_entries(x_bounded)), list())
  val_bounded <- psolve(prob_bounded)
  ## x should be at lower bound (0)
  expect_equal(as.numeric(value(x_bounded)), c(0, 0), tolerance = 1e-3)

  ## Same problem with explicit constraints
  x_constrained <- Variable(2)
  prob_constrained <- Problem(Minimize(sum_entries(x_constrained)),
                              list(x_constrained >= 0, x_constrained <= 5))
  val_constrained <- psolve(prob_constrained)
  expect_equal(as.numeric(value(x_constrained)), c(0, 0), tolerance = 1e-3)

  ## Both approaches give the same optimal value
  expect_equal(val_bounded, val_constrained, tolerance = 1e-3)

  ## Maximize: should hit upper bound
  x_bounded_max <- Variable(2, bounds = list(0, 5))
  prob_max <- Problem(Maximize(sum_entries(x_bounded_max)), list())
  val_max <- psolve(prob_max)
  expect_equal(as.numeric(value(x_bounded_max)), c(5, 5), tolerance = 1e-3)

  ## Bounds with vector lower and upper
  x_vec <- Variable(3, bounds = list(c(1, 2, 3), c(4, 5, 6)))
  prob_vec <- Problem(Minimize(sum_entries(x_vec)), list())
  val_vec <- psolve(prob_vec)
  expect_equal(as.numeric(value(x_vec)), c(1, 2, 3), tolerance = 1e-3)

  ## Bounds with vector lower, scalar upper
  x_mixed <- Variable(2, bounds = list(c(4, 5), 6))
  prob_mixed <- Problem(Maximize(sum_entries(x_mixed)), list())
  val_mixed <- psolve(prob_mixed)
  expect_equal(as.numeric(value(x_mixed)), c(6, 6), tolerance = 1e-3)
})


# =====================================================================
# ADDITIONAL CROSS-CUTTING TESTS
# =====================================================================

# ── NonNeg/NonPos dual parity ────────────────────────────────────────
# CVXPY: test_constraints.py::test_nonneg_dual — dual values match
# between Inequality and NonNeg formulations

## @cvxpy test_constraints.py::TestConstraints::test_nonneg_dual
test_that("NonNeg/Inequality dual value parity (CVXPY parity)", {
  x <- Variable(3)
  c_val <- c(0, 1, 2)

  ## Inequality formulation: c - x <= 0  (same as x >= c)
  prob_ineq <- Problem(Minimize(sum_entries(x)), list(c_val - x <= 0))
  psolve(prob_ineq, solver = CLARABEL_SOLVER)
  dual_ineq <- dual_value(prob_ineq@constraints[[1L]])

  ## NonNeg formulation: x - c >= 0
  x2 <- Variable(3)
  prob_nn <- Problem(Minimize(sum_entries(x2)), list(NonNeg(x2 - c_val)))
  psolve(prob_nn, solver = CLARABEL_SOLVER)
  dual_nn <- dual_value(prob_nn@constraints[[1L]])

  expect_equal(as.numeric(dual_nn), as.numeric(dual_ineq), tolerance = 1e-3)
})

# ── Kron is_dpp behavior ─────────────────────────────────────────────
# CVXPY: kron is NOT DPP when the constant arg is parametric

## @cvxpy NONE
test_that("kron: is_dpp returns FALSE when constant arg has parameters", {
  P <- Parameter(c(2, 2))
  value(P) <- diag(2)
  x <- Variable(c(2, 2))
  expr <- kron(P, x)
  expect_true(is_affine(expr))
  ## kron with Parameter should NOT be DPP (CVXPY limitation)
  expect_false(is_dpp(expr))

  ## kron with Constant IS DPP
  C <- Constant(diag(2))
  expr2 <- kron(C, x)
  expect_true(is_dpp(expr2))
})

# ── Quad form with sparse matrix ─────────────────────────────────────
# CVXPY: test_quad_form.py::test_sparse_quad_form

## @cvxpy test_quad_form.py::TestNonOptimal::test_sparse_quad_form
test_that("quad_form: sparse P matrix", {
  P_sparse <- Matrix::Diagonal(2)
  x <- Variable(2)
  cost <- quad_form(x, P_sparse)
  prob <- Problem(Minimize(cost), list(x == c(1, 2)))
  val <- psolve(prob, solver = "OSQP")
  ## 1^2 + 2^2 = 5
  expect_equal(val, 5, tolerance = 1e-3)
})

# ── Quad form zero matrix ────────────────────────────────────────────
# CVXPY: test_quad_form.py::test_zero_matrix

## @cvxpy test_quad_form.py::TestNonOptimal::test_zero_matrix
test_that("quad_form: zero P reduces to LP", {
  x <- Variable(3)
  P_zero <- matrix(0, 3, 3)
  c_vec <- c(-1, -1, -1)
  ## 0.5 * quad_form(x, 0) + c'x s.t. Ix <= ones
  A <- diag(3)
  b <- rep(1, 3)
  prob <- Problem(Minimize(0.5 * quad_form(x, P_zero) + t(c_vec) %*% x),
                  list(A %*% x <= b))
  val <- psolve(prob, solver = "SCS")
  expect_equal(status(prob), OPTIMAL)
  expect_equal(val, -3.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), c(1, 1, 1), tolerance = 1e-2)
})
