## Solver legs CVXPY runs that CVXR's earlier ports of the same tests dropped.
##
## Found by the 2026-10-07 dropped-solver-leg sweep (plan
## notes/cvxpy_1.9.3_plan.md, Stage 3 finding 2): the upstream test names a
## solver, or several, and the R port of it ran on a different solver or a
## subset. Every leg here was run against the upstream assertion before being
## written down; they guard the legs so they cannot regress unseen. The
## existing ports stay where they are.

## @cvxpy test_conic_solvers.py::TestAllSolvers::test_mixed_integer_behavior
test_that("dropped leg: integer LP on HIGHS", {
  require_solver("HIGHS")
  x <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 0))
  psolve(prob, solver = "HIGHS")
  expect_equal(as.numeric(value(x)), c(0, 0), tolerance = 1e-6)
})

## @cvxpy test_constraints.py::TestConstraints::test_broadcasting
test_that("dropped leg: broadcasting equality on the QP path (OSQP)", {
  require_solver("OSQP")
  x <- Variable(c(3, 1))
  con <- x == Constant(matrix(1, 1L, 3L))
  expect_equal(con@shape, c(3L, 3L))
  psolve(Problem(Minimize(sum_squares(x)), list(con)), solver = "OSQP")
  expect_equal(as.numeric(value(x)), rep(1, 3), tolerance = 1e-4)
})

## @cvxpy test_constraints.py::TestConstraints::test_nonneg_dual
test_that("dropped leg: NonNeg dual on OSQP matches the Inequality dual", {
  require_solver("OSQP"); require_solver("CLARABEL")
  cc <- 0:2
  x <- Variable(3)
  p0 <- Problem(Minimize(sum_entries(x)), list(cc - x <= 0))
  psolve(p0, solver = "CLARABEL")
  ref <- dual_value(p0@constraints[[1L]])
  x1 <- Variable(3)
  p1 <- Problem(Minimize(sum_entries(x1)), list(NonNeg(x1 - cc)))
  psolve(p1, solver = "OSQP")
  expect_equal(as.numeric(dual_value(p1@constraints[[1L]])), as.numeric(ref),
               tolerance = 1e-4)
})

## TestMIPVariable bool_prob / int_prob / bool_socp / int_socp.
.dl_mip_leg <- function(solver) {
  xb <- Variable(boolean = TRUE)
  expect_equal(psolve(Problem(Minimize(abs(xb - 0.2))), solver = solver), 0.2,
               tolerance = 1e-5, info = solver)
  expect_equal(as.numeric(value(xb)), 0, tolerance = 1e-5, info = solver)
  t <- Variable(); xb <- Variable(boolean = TRUE)
  expect_equal(psolve(Problem(Minimize(t), list(abs(xb) <= t)), solver = solver), 0,
               tolerance = 1e-5, info = solver)
  C <- matrix(c(0, 1, 0, 1, 1, 1), 3L, 2L)
  A <- Variable(c(3, 2), boolean = TRUE)
  expect_equal(psolve(Problem(Minimize(sum_entries(abs(A - C)))), solver = solver), 0,
               tolerance = 1e-5, info = solver)
  expect_equal(value(A), C, tolerance = 1e-4, info = solver)
  t <- Variable(); A <- Variable(c(3, 2), boolean = TRUE)
  expect_equal(psolve(Problem(Minimize(t), list(sum_entries(abs(A - C)) <= t)),
                      solver = solver), 0, tolerance = 1e-5, info = solver)
  yi <- Variable(integer = TRUE)
  expect_equal(psolve(Problem(Minimize(abs(yi - 0.2))), solver = solver), 0.2,
               tolerance = 1e-5, info = solver)
  yi <- Variable(integer = TRUE); t <- Variable()
  p <- Problem(Minimize(t), list(yi == 0.5, t >= 0))
  ## GUROBI reports INF_OR_UNBD here and CVXR warns so; upstream accepts that.
  withCallingHandlers(psolve(p, solver = solver), warning = function(w) {
    if (grepl("infeasible or unbounded", conditionMessage(w))) invokeRestart("muffleWarning")
  })
  expect_true(status(p) %in% c("infeasible", "unbounded", "infeasible_or_unbounded",
                               "infeasible_inaccurate", "unbounded_inaccurate"),
              info = solver)
  if (solver != "XPRESS") {   # upstream issue #1815
    xb <- Variable(boolean = TRUE); t <- Variable()
    expect_equal(psolve(Problem(Minimize(t), list(square(xb - 0.2) <= t)), solver = solver),
                 0.04, tolerance = 1e-4, info = solver)
  }
  yi <- Variable(integer = TRUE); t <- Variable()
  expect_equal(psolve(Problem(Minimize(t), list(square(yi - 0.2) <= t)), solver = solver),
               0.04, tolerance = 1e-4, info = solver)
  expect_equal(as.numeric(value(yi)), 0, tolerance = 1e-5, info = solver)
}

## @cvxpy test_mip_vars.py::TestMIPVariable::test_all_solvers
test_that("dropped leg: MIP variables on CPLEX, GUROBI, MOSEK, XPRESS", {
  ran <- 0L
  for (solver in c("CPLEX", "GUROBI", "MOSEK", "XPRESS")) {
    if (!solver %in% installed_solvers()) next
    .dl_mip_leg(solver)
    ran <- ran + 1L
  }
  if (ran == 0L) skip("none of CPLEX, GUROBI, MOSEK, XPRESS installed")
})

## @cvxpy test_pow_cone_nd.py::TestPowConeND::test_3d_pow_cone_scalar_alpha
test_that("dropped leg: PowCone3D with scalar alpha on SCS", {
  require_solver("SCS")
  x <- Variable(3)
  prob <- Problem(Minimize(p_norm(x, 2)), list(PowCone3D(x[1], x[2], x[3], 0.75)))
  psolve(prob, solver = "SCS")
  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
})

## TestProblem.test_dual_variables. CVXPY's nested list [[1, 2], [3, 4]] is
## read as COLUMNS there, so the matrix is matrix(1:4, 2) in R.
.dl_dual_leg <- function(solver, tol) {
  x <- Variable(2); z <- Variable(2)
  p <- Problem(Minimize(norm1(x + z)),
               list(x >= c(2, 3), matrix(c(1, 2, 3, 4), 2L) %*% z == c(-1, -4),
                    p_norm(x + z, 2) <= 100))
  expect_equal(psolve(p, solver = solver), 4, tolerance = tol, info = solver)
  expect_equal(as.numeric(value(x)), c(4, 3), tolerance = tol, info = solver)
  expect_equal(as.numeric(value(z)), c(-4, 1), tolerance = tol, info = solver)
  expect_lt(max(abs(as.numeric(dual_value(p@constraints[[1L]])) - c(0, 1))), tol)
  expect_lt(max(abs(as.numeric(dual_value(p@constraints[[2L]])) - c(-1, 0.5))), tol)
  expect_lt(abs(as.numeric(dual_value(p@constraints[[3L]]))), tol)

  A <- Variable(c(2, 2)); B <- Variable(c(2, 2)); C <- Variable(c(3, 2))
  TT <- matrix(2, 2L, 3L)
  p <- Problem(Minimize(1), list(A >= TT %*% C, A == B, C == t(TT)))
  psolve(p, solver = solver)
  expect_lt(max(abs(as.numeric(dual_value(p@constraints[[1L]])))), tol)
  expect_lt(max(abs(as.numeric(dual_value(p@constraints[[2L]])))), tol)
  expect_lt(max(abs(as.numeric(dual_value(p@constraints[[3L]])))), tol)
}

## @cvxpy test_problem.py::TestProblem::test_dual_variables
test_that("dropped leg: dual variable recovery on SCS and CVXOPT", {
  require_solver("SCS"); require_solver("CVXOPT")
  .dl_dual_leg("SCS", 0.05)       # upstream: 1 decimal place for SCS
  .dl_dual_leg("CVXOPT", 1e-5)
})

## @cvxpy test_problem.py::TestProblem::test_get_problem_data
test_that("dropped leg: problem data for CVXOPT", {
  require_solver("CVXOPT")
  x <- Variable(2)
  data <- problem_data(Problem(Minimize(p_norm(x, 2) + 3)), solver = "CVXOPT")$data
  expect_equal(data$dims@soc, 3L)
})

## @cvxpy test_problem.py::TestProblem::test_invalid_solvers
test_that("dropped leg: SCS refuses an unbounded objective it cannot solve", {
  require_solver("SCS")
  a <- Variable()
  expect_error(psolve(Problem(Minimize(a)), solver = "SCS"))
})

## @cvxpy test_problem.py::TestProblem::test_psd_duals
test_that("dropped leg: PSD duals on CVXOPT", {
  require_solver("CVXOPT")
  C <- Variable(c(2, 2), symmetric = TRUE)
  cons <- list(C %<<% diag(2, 2L))
  expect_equal(psolve(Problem(Maximize(C[1, 1]), cons), solver = "CVXOPT"), 2,
               tolerance = 1e-6)
  ref <- dual_value(cons[[1L]])
  C <- Variable(c(2, 2), symmetric = TRUE)
  X <- Variable(c(2, 2), PSD = TRUE)
  cons2 <- list(X == diag(2, 2L) - C)
  psolve(Problem(Maximize(C[1, 1]), cons2), solver = "CVXOPT")
  d <- dual_value(cons2[[1L]])
  expect_equal((d + t(d)) / 2, ref, tolerance = 1e-5)
})

## @cvxpy test_problem.py::TestProblem::test_scalar_lp
test_that("dropped leg: unbounded scalar LP on CVXOPT", {
  require_solver("CVXOPT")
  a <- Variable()
  p <- Problem(Minimize(-a), list(a >= 2))
  r <- psolve(p, solver = "CVXOPT")
  expect_equal(status(p), "unbounded")
  expect_true(is.infinite(r) && r < 0)
})

## @cvxpy test_problem.py::TestProblem::test_solver_stats
test_that("dropped leg: SCS solver stats, including extra_stats$info", {
  require_solver("SCS")
  x <- Variable(2)
  p <- Problem(Minimize(p_norm(x, 2)), list(x == 0))
  psolve(p, solver = "SCS")
  st <- solver_stats(p)
  expect_gt(st@solve_time, 0)
  expect_gt(st@setup_time, 0)
  expect_gt(st@num_iters, 0)
  expect_false(is.null(st@extra_stats[["info"]]))
})

## @cvxpy test_qp_solvers.py::TestQp::test_highs_dense_quad_form
test_that("dropped leg: dense quad_form on HIGHS agrees with CLARABEL", {
  require_solver("HIGHS"); require_solver("CLARABEL")
  set.seed(42)
  nv <- 60L; nn <- 20L
  M <- matrix(0, nv, nn)
  for (i in seq_len(nv)) M[i, sample(nn, 2L)] <- c(1, -1)
  B <- matrix(rnorm(nn * nn), nn)
  Sigma <- tcrossprod(B)              # dense PSD
  vals <- vapply(c("HIGHS", "CLARABEL"), function(solver) {
    x <- Variable(nv, nonneg = TRUE)
    y <- t(M) %*% x
    psolve(Problem(Maximize(sum_entries(x) - 0.1 * quad_form(y, Sigma)), list(x <= 10)),
           solver = solver)
  }, numeric(1))
  expect_equal(vals[["HIGHS"]], vals[["CLARABEL"]], tolerance = 1e-4)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_entr
test_that("dropped leg: entr on SCS", {
  require_solver("SCS")
  for (n in c(5L, 10L, 25L)) {
    x <- Variable(n)
    psolve(Problem(Maximize(sum(entr(x))), list(sum(x) == 1)), solver = "SCS")
    expect_equal(as.numeric(value(x)), rep(1 / n, n), tolerance = 1e-4, info = n)
  }
})

## @cvxpy test_conic_solvers.py::TestSCS::test_log_problem
test_that("dropped leg: log in objective, constraint, and indexed, on SCS", {
  require_solver("SCS")
  x <- Variable(2)
  expect_equal(psolve(Problem(Maximize(sum(log(x))), list(x <= c(1, exp(1)))),
                      solver = "SCS"), 1, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, exp(1)), tolerance = 1e-4)
  x <- Variable(2)
  expect_equal(psolve(Problem(Minimize(sum(x)), list(log(x) >= 0, x <= c(1, 1))),
                      solver = "SCS"), 2, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-4)
  x <- Variable(2)
  expect_equal(psolve(Problem(Maximize(log(x)[2]), list(x <= c(1, exp(1)))),
                      solver = "SCS"), 1, tolerance = 1e-4)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_solve_problem_twice
test_that("dropped leg: solving twice on SCS gives the same answer", {
  require_solver("SCS")
  x <- Variable(5)
  p <- Problem(Maximize(sum(log(x))), list(sum(x) == 1))
  psolve(p, solver = "SCS")
  first <- as.numeric(value(x))
  expect_equal(first, rep(0.2, 5), tolerance = 1e-4)
  psolve(p, solver = "SCS")
  expect_equal(as.numeric(value(x)), first, tolerance = 1e-6)
})

## @cvxpy test_bounds.py::TestHIGHSBoundsPropagation::test_get_problem_data_array_bounds
test_that("dropped leg: per-element bounds reach HIGHS problem data", {
  require_solver("HIGHS")
  x <- Variable(3, bounds = list(c(1, 2, 3), c(10, 20, 30)))
  data <- problem_data(Problem(Minimize(sum_entries(x))), solver = "HIGHS")$data
  expect_equal(as.numeric(data[[LOWER_BOUNDS]]), c(1, 2, 3))
  expect_equal(as.numeric(data[[UPPER_BOUNDS]]), c(10, 20, 30))
})
