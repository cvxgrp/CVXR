## CVXPY 1.9 infeasibility-certificate parity
## (solver_test_helpers.py::StandardTestInfeasibleProblems).
##
## CVXPY v1.9.0 (#3197 + #3228) propagates a solver's certificate of
## infeasibility (the dual ray) up to each constraint's dual_value, even though
## the problem has no primal solution. CVXR mirrors this for CLARABEL: the
## certificate flows solver-invert -> failure_solution -> ConeMatrixStuffing
## invert (hoisted above the failure short-circuit) -> problem unpack.
##
## verify_post_solve_guarantees mirrors the CVXPY helper: infeasible status,
## +/-Inf value, NULL variable values, and a non-NULL dual per constraint
## (every element for composite SOC/ExpCone/PowCone duals).

.infeasible_guarantees <- function(prob, vars, cons, maximize = FALSE,
                                   solver = "CLARABEL", check_duals = TRUE) {
  psolve(prob, solver = solver)
  expect_equal(status(prob), "infeasible")
  expect_equal(as.numeric(value(prob)), if (maximize) -Inf else Inf)
  for (v in vars) expect_null(value(v))
  if (!check_duals) return(invisible())
  for (con in cons) {
    dv <- dual_value(con)
    expect_false(is.null(dv))
    parts <- if (is.list(dv)) dv else list(dv)
    for (p in parts) expect_false(is.null(p))
  }
}

## @cvxpy solver_test_helpers.py::StandardTestInfeasibleProblems::test_lp_ineq_constraints
test_that("infeasible LP (inequalities) yields a Farkas certificate dual", {
  A <- matrix(c(1, 0, 0, 1, -1, -1), nrow = 3, byrow = TRUE)
  b <- c(0, 0, -1)
  x <- Variable(2)
  con <- A %*% x <= b
  prob <- Problem(Minimize(0), list(con))
  .infeasible_guarantees(prob, list(x), list(con))

  ## Farkas certificate: y >= 0, A^T y == 0, b^T y < 0.
  y <- as.numeric(dual_value(con))
  expect_true(all(y >= -1e-6))
  expect_equal(as.numeric(t(A) %*% y), c(0, 0), tolerance = 1e-4)
  expect_true(sum(b * y) < 0)
})

## @cvxpy solver_test_helpers.py::StandardTestInfeasibleProblems::test_lp_eq_constraints
test_that("infeasible LP (equalities) yields a Farkas certificate dual", {
  A <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, byrow = TRUE)
  b <- c(0, 0, 1)
  x <- Variable(2)
  ceq <- A %*% x == b
  cnn <- x >= 0
  prob <- Problem(Minimize(0), list(ceq, cnn))
  .infeasible_guarantees(prob, list(x), list(ceq, cnn))

  ## Farkas certificate: A^T y >= 0, b^T y < 0.
  y <- as.numeric(dual_value(ceq))
  expect_true(min(as.numeric(t(A) %*% y)) >= -1e-4)
  expect_true(sum(b * y) < 0)
})

## @cvxpy solver_test_helpers.py::StandardTestInfeasibleProblems::test_soc
test_that("infeasible SOC yields a certificate dual", {
  x <- Variable(2)
  con <- SOC(-1, x)
  prob <- Problem(Minimize(0), list(con))
  .infeasible_guarantees(prob, list(x), list(con))
})

## @cvxpy solver_test_helpers.py::StandardTestInfeasibleProblems::test_power_cone_3d
test_that("infeasible PowCone3D yields a certificate dual", {
  z <- Variable()
  con <- PowCone3D(-1, 1, z, 0.5)
  prob <- Problem(Minimize(0), list(con))
  .infeasible_guarantees(prob, list(z), list(con))
})

## @cvxpy solver_test_helpers.py::StandardTestInfeasibleProblems::test_power_cone_nd
test_that("infeasible PowConeND yields a certificate dual", {
  x <- Variable(3)
  z <- Variable()
  cpn <- PowConeND(x, z, c(0.2, 0.3, 0.5))
  c2 <- x == rep(1, 3)
  c3 <- z == 2.0
  prob <- Problem(Minimize(0), list(cpn, c2, c3))
  .infeasible_guarantees(prob, list(x, z), list(cpn, c2, c3))
})

## @cvxpy solver_test_helpers.py::StandardTestInfeasibleProblems::test_exp_cone
test_that("infeasible ExpCone yields a certificate dual", {
  x <- Variable(); y <- Variable(); z <- Variable()
  ce <- ExpCone(x, y, z)
  c2 <- x == 1; c3 <- y == 1; c4 <- z == 1
  prob <- Problem(Minimize(0), list(ce, c2, c3, c4))
  .infeasible_guarantees(prob, list(x, y, z), list(ce, c2, c3, c4))
})

## @cvxpy solver_test_helpers.py::StandardTestInfeasibleProblems::test_psd_cone
test_that("infeasible PSD cone yields a certificate dual", {
  X <- Variable(c(2, 2))
  c1 <- PSD(X)
  c2 <- X[1, 1] == -1
  prob <- Problem(Minimize(0), list(c1, c2))
  .infeasible_guarantees(prob, list(X), list(c1, c2))
})

## @cvxpy solver_test_helpers.py::StandardTestInfeasibleProblems::test_soc_exp_mixed
test_that("infeasible SOC + ExpCone mixed yields certificate duals", {
  x <- Variable(); y <- Variable()
  c1 <- x == y
  c2 <- SOC(x, hstack(1))       # x >= 1
  c3 <- ExpCone(x, y, 1)        # x <= 0
  prob <- Problem(Maximize(0), list(c1, c2, c3))
  .infeasible_guarantees(prob, list(x, y), list(c1, c2, c3), maximize = TRUE)
})

## ---------------------------------------------------------------------
## QP-path solvers: test_qp_solvers.py::TestQp wrappers around the same
## StandardTestInfeasibleProblems LP helpers, parametrized by solver.
## ---------------------------------------------------------------------

## OSQP exposes prim_inf_cert (the primal-infeasibility certificate) on
## infeasible solves; CVXR v1.9.0 #3228 plumbs it through OSQP invert ->
## failure_solution -> problem unpack. Full Farkas check.

## @cvxpy test_qp_solvers.py::TestQp::test_osqp_infeasible_lp_ineq_constraints
test_that("OSQP infeasible LP (inequalities) yields a Farkas certificate dual", {
  skip_if_not_installed("osqp")
  A <- matrix(c(1, 0, 0, 1, -1, -1), nrow = 3, byrow = TRUE)
  b <- c(0, 0, -1)
  x <- Variable(2)
  con <- A %*% x <= b
  prob <- Problem(Minimize(0), list(con))
  .infeasible_guarantees(prob, list(x), list(con), solver = "OSQP")

  ## Farkas certificate: y >= 0, A^T y == 0, b^T y < 0.
  y <- as.numeric(dual_value(con))
  expect_true(all(y >= -1e-6))
  expect_equal(as.numeric(t(A) %*% y), c(0, 0), tolerance = 1e-3)
  expect_true(sum(b * y) < 0)
})

## @cvxpy test_qp_solvers.py::TestQp::test_osqp_infeasible_lp_eq_constraints
test_that("OSQP infeasible LP (equalities) yields a Farkas certificate dual", {
  skip_if_not_installed("osqp")
  A <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, byrow = TRUE)
  b <- c(0, 0, 1)
  x <- Variable(2)
  ceq <- A %*% x == b
  cnn <- x >= 0
  prob <- Problem(Minimize(0), list(ceq, cnn))
  .infeasible_guarantees(prob, list(x), list(ceq, cnn), solver = "OSQP")

  ## Farkas certificate: A^T y >= 0, b^T y < 0.
  y <- as.numeric(dual_value(ceq))
  expect_true(min(as.numeric(t(A) %*% y)) >= -1e-4)
  expect_true(sum(b * y) < 0)
})

## HiGHS reports infeasible status / +Inf value / NULL primal correctly,
## but the R `highs` package (1.12) exposes no getDualRay(), so the
## infeasibility certificate (CVXPY's highs_conif.py #3228) cannot be
## recovered. We verify the satisfiable guarantees, then skip the
## certificate assertion as @v19-pending until R `highs` adds a dual-ray
## accessor (same upstream-capability gap as HiGHS warm-start).

## @cvxpy test_qp_solvers.py::TestQp::test_highs_infeasible_lp_ineq_constraints
## @v19-pending: highs-infeasibility-certificate [now]
test_that("HiGHS infeasible LP (inequalities): status/value/primal (cert pending)", {
  skip_if_not_installed("highs")
  A <- matrix(c(1, 0, 0, 1, -1, -1), nrow = 3, byrow = TRUE)
  b <- c(0, 0, -1)
  x <- Variable(2)
  con <- A %*% x <= b
  prob <- Problem(Minimize(0), list(con))
  .infeasible_guarantees(prob, list(x), list(con),
                         solver = "HIGHS", check_duals = FALSE)
  skip("@v19-pending: HiGHS infeasibility certificate needs R highs getDualRay() (unimplemented upstream)")
})

## @cvxpy test_qp_solvers.py::TestQp::test_highs_infeasible_lp_eq_constraints
## @v19-pending: highs-infeasibility-certificate [now]
test_that("HiGHS infeasible LP (equalities): status/value/primal (cert pending)", {
  skip_if_not_installed("highs")
  A <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, byrow = TRUE)
  b <- c(0, 0, 1)
  x <- Variable(2)
  ceq <- A %*% x == b
  cnn <- x >= 0
  prob <- Problem(Minimize(0), list(ceq, cnn))
  .infeasible_guarantees(prob, list(x), list(ceq, cnn),
                         solver = "HIGHS", check_duals = FALSE)
  skip("@v19-pending: HiGHS infeasibility certificate needs R highs getDualRay() (unimplemented upstream)")
})

## @cvxpy test_qp_solvers.py::TestQp::test_highs_dense_quad_form
test_that("HiGHS solves a dense quad_form on a linear expression (issue #3301)", {
  skip_if_not_installed("highs")
  ## CVXPY #3301/#3302: a dense quad_form on a linear expression yields a
  ## Hessian whose two triangles differ by fp-epsilon. HiGHS >= 1.14 adds a
  ## strict symmetry check and REJECTS the full-Q form (aborts in
  ## hi_new_solver() with "Square Hessian contains N non-symmetries";
  ## reproduced against the 1.14 dual-ray fork). HiGHS 1.12 -- the current
  ## CRAN / supported target -- has no such check. highs_qpif.R now
  ## symmetrizes Q (`(Q + t(Q))/2`, the exact canonical Hessian) so the model
  ## is accepted on BOTH 1.12 and 1.14. (CVXPY instead passes a triangle via
  ## sp.triu + kTriangular; CVXR cannot, because the R highs API exposes only
  ## the "square" Hessian format -- see the deficiency note in highs_qpif.R.)
  ## Smoke test: the problem must solve to OPTIMAL.
  set.seed(42)
  n_vars <- 60L; n_nodes <- 20L
  rows <- integer(0); cols <- integer(0); vals <- numeric(0)
  for (i in seq_len(n_vars)) {
    jk <- sample(n_nodes, 2)
    rows <- c(rows, i, i); cols <- c(cols, jk); vals <- c(vals, 1, -1)
  }
  M <- as.matrix(Matrix::sparseMatrix(i = rows, j = cols, x = vals,
                                      dims = c(n_vars, n_nodes)))
  A   <- matrix(rnorm(n_nodes * n_nodes), n_nodes)
  raw <- A %*% t(A) + matrix(rnorm(n_nodes * n_nodes), n_nodes) * 0.01
  sym <- 0.5 * (raw + t(raw))
  e   <- eigen(sym, symmetric = TRUE)
  Sigma <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors)

  x <- Variable(n_vars, nonneg = TRUE)
  y <- t(M) %*% x                       # x @ M  ->  (n_nodes, 1)
  prob <- Problem(Maximize(sum(x) - 0.1 * quad_form(y, Sigma)), list(x <= 10))
  psolve(prob, solver = "HIGHS")
  expect_equal(status(prob), "optimal")
})
