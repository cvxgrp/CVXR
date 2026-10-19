## Claims audit, decision 24c (CVXR 1.9.2.9160): the remaining solver rows.
##
##   C0742  base `reduction_invert(ConicSolver)` adds the objective OFFSET
##          (conic_solver.py:321); latent, every shipped conic solver overrides.
##   C0900  the QP guard rejects ND power cones too (qp_solver.py:28-31).
##   C0886  EXTRA_STATS on every solver CVXPY 1.9.3 sets it for and CVXR did
##          not: OSQP (osqp_qpif.py:66), ECOS (ecos_conif.py:164), ECOS_BB
##          (ecos_bb_conif.py:81), Gurobi conic and QP (gurobi_conif.py:114,
##          gurobi_qpif.py:99), DIFFCP (diffcp_conif.py:102); OSQP's NUM_ITERS
##          moves inside the solution branch (osqp_qpif.py:78).
##   C0563  `reduction_invert(Canonicalization)` filters primal_vars to the
##          variables of the problem it received (canonicalization.py:78-79).
##   C0892  PIQP `backend`: validated (piqp_qpif.py:92-93), a dense request
##          builds a DENSE model (`.toarray()`, :164-166; on an in-place update
##          :120/:126/:132), and a cached model of the other backend is rebuilt (:144-147).
##   C0730  `accept_unknown`: CVXR promotes InsufficientProgress only when the
##          option is TRUE; CVXPY on key presence. Kept (ADR D_PORT.12); the
##          carrier is pinned here. A pin of existing behavior: it cannot be
##          seen failing.
##
## Every other test in this file was seen failing on 1.9.2.9159 first.

# -- C0742: base ConicSolver invert ---------------------------------------

## @cvxpy test_conic_solvers.py::test_offset_in_opt_val
test_that("base reduction_invert(ConicSolver) adds the objective offset", {
  slv <- CVXR:::ConicSolver()
  sol <- list(status = "optimal", value = 2.5, primal = c(1, 1.5),
              eq_dual = numeric(0), ineq_dual = numeric(0))
  inv <- list(var_id = 7L, offset = 1000, eq_constr = list(),
              other_constr = list())
  out <- CVXR:::reduction_invert(slv, sol, inv)
  expect_equal(out@status, "optimal")
  expect_equal(out@opt_val, 1002.5)
  expect_equal(out@primal_vars[["7"]], c(1, 1.5))
  ## A failure carries no value to offset.
  bad <- CVXR:::reduction_invert(slv, list(status = "infeasible"), inv)
  expect_equal(bad@status, "infeasible")
})

# -- C0900: QP guard, ND power cone ---------------------------------------

## @cvxpy test_qp_solvers.py::TestQpSolverValidation::test_qp_solver_rejects_soc_cones
test_that("QpSolver apply rejects ND power cones like the other cones", {
  mk_dims <- function(prop, value) {
    dims <- CVXR:::ConeDims(list())
    S7::prop(dims, prop) <- value
    dims
  }
  slv <- CVXR:::OSQP_QP_Solver()
  expect_error(
    CVXR:::reduction_apply(slv, list(x_id = 1L, dims = mk_dims("pnd", list(c(0.3, 0.7))))),
    "cannot handle"
  )
  ## Existing guards, unchanged.
  expect_error(
    CVXR:::reduction_apply(slv, list(x_id = 1L, dims = mk_dims("p3d", 0.5))),
    "cannot handle"
  )
  expect_error(
    CVXR:::reduction_apply(slv, list(x_id = 1L, dims = mk_dims("soc", 3L))),
    "cannot handle"
  )
})

# -- C0886: EXTRA_STATS, the systemic set ---------------------------------

## @cvxpy NONE
test_that("OSQP: solver_stats carries the solver result as extra_stats", {
  sth <- sth_qp_0()
  psolve(sth$prob, solver = "OSQP")
  ss <- solver_stats(sth$prob)
  expect_false(is.null(ss@extra_stats))
  expect_false(is.null(ss@extra_stats$info$status_val))
  expect_true(is.null(ss@extra_stats$.len_eq))
  expect_false(is.null(ss@num_iters))
})

## @cvxpy NONE
test_that("OSQP: num_iters is reported only with a solution, as CVXPY", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x)), list(x >= 1, x <= 0))
  psolve(prob, solver = "OSQP")
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
  ss <- solver_stats(prob)
  expect_true(is.null(ss@num_iters))
  expect_false(is.null(ss@extra_stats))
})

## @cvxpy NONE
test_that("ECOS: solver_stats carries the solver result as extra_stats", {
  require_solver("ECOS")
  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x - c(1, 2), 2)), list(sum(x) <= 1))
  psolve(prob, solver = "ECOS")
  ss <- solver_stats(prob)
  expect_false(is.null(ss@extra_stats))
  expect_false(is.null(ss@extra_stats$retcodes))
})

## @cvxpy NONE
test_that("ECOS_BB: solver_stats carries the solver result as extra_stats", {
  require_solver("ECOS_BB")
  sth <- sth_mi_lp_0()
  psolve(sth$prob, solver = "ECOS_BB")
  ss <- solver_stats(sth$prob)
  expect_false(is.null(ss@extra_stats))
  expect_false(is.null(ss@extra_stats$retcodes))
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_warm_start
test_that("Gurobi conic and QP: solver_stats carries the Gurobi result as extra_stats", {
  require_solver("GUROBI")
  sth <- sth_mi_lp_0()
  psolve(sth$prob, solver = "GUROBI")
  ss <- solver_stats(sth$prob)
  expect_false(is.null(ss@extra_stats))
  expect_equal(ss@extra_stats$status, "OPTIMAL")
  expect_true(is.null(ss@extra_stats$.n_orig))

  sth <- sth_qp_0()
  psolve(sth$prob, solver = "GUROBI")
  ss <- solver_stats(sth$prob)
  expect_false(is.null(ss@extra_stats))
  expect_equal(ss@extra_stats$status, "OPTIMAL")
  expect_true(is.null(ss@extra_stats$.len_eq))
})

## @cvxpy NONE
test_that("DIFFCP: solver_stats carries the diffcp result as extra_stats", {
  skip_if_not_installed("diffcp")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - c(1, 2))), list(sum(x) <= 1))
  psolve(prob, requires_grad = TRUE)
  ss <- solver_stats(prob)
  expect_false(is.null(ss@extra_stats))
  expect_false(is.null(ss@extra_stats$info))
})

# -- C0563: Canonicalization invert filters primal_vars -------------------

## @cvxpy NONE
test_that("solution(prob) holds only the problem's own variables", {
  x <- Variable(3)
  prob <- Problem(Minimize(norm1(x - c(1, 2, 3)) + sum_squares(x)), list(x >= 0.5))
  psolve(prob, solver = "CLARABEL")
  sol <- solution(prob)
  expect_equal(names(sol@primal_vars), as.character(x@id))
  expect_equal(as.numeric(sol@primal_vars[[1L]]), as.numeric(value(x)))

  ## Unconstrained: no original-constraint duals to remap, same filter.
  y <- Variable(3)
  prob2 <- Problem(Minimize(norm1(y - c(1, 2, 3)) + sum_squares(y)))
  psolve(prob2, solver = "CLARABEL")
  expect_equal(names(solution(prob2)@primal_vars), as.character(y@id))
})

# -- C0892: PIQP backend ---------------------------------------------------

## @cvxpy test_qp_solvers.py::TestPiqpInterface::test_piqp_dense_backend
test_that("PIQP: the dense backend agrees with sparse and builds a dense model", {
  require_solver("PIQP")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - 3)), list(x[1] + x[2] <= 4))
  psolve(prob, solver = "PIQP", backend = "dense")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), c(2, 2), tolerance = 1e-4)
  model <- get("PIQP", envir = CVXR:::.get_solver_cache(prob))$model
  expect_true(model@dense_backend)
})

## @cvxpy test_qp_solvers.py::TestPiqpInterface::test_piqp_dense_backend
test_that("PIQP: a backend switch under warm_start rebuilds the model", {
  require_solver("PIQP")
  cached_dense <- function(prob) {
    get("PIQP", envir = CVXR:::.get_solver_cache(prob))$model@dense_backend
  }
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - 3)), list(x[1] + x[2] <= 4))
  psolve(prob, solver = "PIQP")
  expect_false(cached_dense(prob))
  psolve(prob, solver = "PIQP", backend = "dense", warm_start = TRUE)
  expect_true(cached_dense(prob))
  expect_equal(as.numeric(value(x)), c(2, 2), tolerance = 1e-4)
  ## Same backend again: the dense model is updated in place, not rebuilt.
  psolve(prob, solver = "PIQP", backend = "dense", warm_start = TRUE)
  expect_true(cached_dense(prob))
  psolve(prob, solver = "PIQP", backend = "sparse", warm_start = TRUE)
  expect_false(cached_dense(prob))
  expect_equal(as.numeric(value(x)), c(2, 2), tolerance = 1e-4)
})

# -- C0730: accept_unknown carrier, pinned --------------------------------

## @cvxpy NONE
test_that("CLARABEL: accept_unknown = FALSE is carried as FALSE (recorded divergence)", {
  x <- Variable(nonneg = TRUE)
  prob <- Problem(Minimize(sum(x)), list(x >= 3))
  pd <- problem_data(prob, solver = "CLARABEL")
  raw <- solve_via_data(pd$chain, pd$data, warm_start = FALSE, verbose = FALSE,
                        solver_opts = list(accept_unknown = FALSE))
  expect_false(raw$accept_unknown)
  raw <- solve_via_data(pd$chain, pd$data, warm_start = FALSE, verbose = FALSE,
                        solver_opts = list(accept_unknown = TRUE))
  expect_true(raw$accept_unknown)
})

# -- C0739: ConicSolver accepts refuses a MIP for a non-MIP solver ------------
## conic_solver.py:121 (MIP_CAPABLE or not problem.is_mixed_integer()). Measured
## on CVXPY 1.9.3: Chain.accepts and CLARABEL.accepts are False for this MILP
## (notes/audit/repro/29). SCIP and ECOS_BB use the base method and are
## MIP-capable; accepts does not call the solver, so no package is needed.

## @cvxpy NONE
test_that("ConicSolver accepts: a MILP is refused unless the solver is MIP-capable", {
  lp <- function(integer) {
    x <- Variable(2, integer = integer)
    Problem(Minimize(x[1] + x[2]), list(x >= 0.5, x <= 3))
  }
  front <- list(Dcp2Cone(), CvxAttr2Constr(), ConeMatrixStuffing())
  for (name in c("CLARABEL", "SCS", "SCIP", "ECOS_BB")) {
    solver <- SOLVER_MAP_CONIC[[name]]
    for (integer in c(FALSE, TRUE)) {
      prob <- lp(integer)
      stuffed <- reduction_apply(Chain(reductions = front), prob)[[1L]]
      expected <- !integer || solver@MIP_CAPABLE
      expect_identical(reduction_accepts(solver, stuffed), expected,
                       info = paste(name, "integer =", integer))
      expect_identical(reduction_accepts(Chain(reductions = c(front, list(solver))), prob),
                       expected, info = paste(name, "chain, integer =", integer))
    }
  }
})
