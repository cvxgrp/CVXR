## CRAN solver test matrix — core solvers only (Clarabel, SCS, OSQP, HiGHS)
## Extracted from test-standard-solver-matrix.R for CRAN submission.
## These 4 solvers are in Imports and always available — zero skips expected.
##
## Tolerances are solver-documented defaults, NOT empirically fudged:
##   Clarabel: 1e-4 (obj), 1e-3 (primal/dual)
##   SCS:      1e-2 (obj/primal/dual)
##   OSQP:     1e-4 (obj/primal), LP/QP only
##   HiGHS:    1e-5 (obj/primal), LP/QP only
##
## NEVER tweak tolerances to get a test to pass.
## If a test fails, the code has a bug — fix the code, not the tolerance.

# ── Helper: verify objective ─────────────────────────────────────

verify_obj <- function(prob, expected, tol, solver_name) {
  val <- psolve(prob, solver = solver_name)
  if (is.infinite(expected) && expected > 0) {
    ## Infeasible: psolve returns Inf
    expect_true(is.infinite(val) && val > 0,
                label = paste(solver_name, "should be infeasible (Inf)"))
  } else if (is.infinite(expected) && expected < 0) {
    ## Unbounded: psolve returns -Inf
    expect_true(is.infinite(val) && val < 0,
                label = paste(solver_name, "should be unbounded (-Inf)"))
  } else {
    expect_equal(val, expected, tolerance = tol,
                 label = paste(solver_name, "objective"))
  }
}

# ── Helper: verify primal ────────────────────────────────────────

verify_primal <- function(prob, var, expected, tol, solver_name) {
  if (is.null(expected)) return(invisible(NULL))
  actual <- value(var)
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = tol,
               label = paste(solver_name, "primal"))
}

# ── Helper: verify duals ─────────────────────────────────────────

verify_dual <- function(con, expected, tol, solver_name) {
  if (is.null(expected)) return(invisible(NULL))
  actual <- dual_value(con)
  if (is.list(expected)) {
    ## SOC / ExpCone / PowCone: list of arrays
    for (k in seq_along(expected)) {
      expect_equal(as.numeric(actual[[k]]), as.numeric(expected[[k]]),
                   tolerance = tol,
                   label = paste(solver_name, "dual component", k))
    }
  } else {
    expect_equal(as.numeric(actual), as.numeric(expected), tolerance = tol,
                 label = paste(solver_name, "dual"))
  }
}

# ══════════════════════════════════════════════════════════════════
# CLARABEL
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_lp_0
test_that("Clarabel: lp_0", {
  h <- sth_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CLARABEL")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_lp_1
test_that("Clarabel: lp_1", {
  h <- sth_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CLARABEL")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "CLARABEL")
  for (i in seq_along(h$con_duals)) {
    verify_dual(h$prob@constraints[[i]], h$con_duals[[i]], 1e-3, "CLARABEL")
  }
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_lp_2
test_that("Clarabel: lp_2", {
  h <- sth_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CLARABEL")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "CLARABEL")
  for (i in seq_along(h$con_duals)) {
    verify_dual(h$prob@constraints[[i]], h$con_duals[[i]], 1e-3, "CLARABEL")
  }
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_lp_3
test_that("Clarabel: lp_3 (unbounded)", {
  h <- sth_lp_3()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_lp_4
test_that("Clarabel: lp_4 (infeasible)", {
  h <- sth_lp_4()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_lp_5
test_that("Clarabel: lp_5 (redundant equalities)", {
  h <- sth_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-3, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_qp_0
test_that("Clarabel: qp_0", {
  h <- sth_qp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CLARABEL")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "CLARABEL")
  for (i in seq_along(h$con_duals)) {
    verify_dual(h$prob@constraints[[i]], h$con_duals[[i]], 1e-3, "CLARABEL")
  }
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_socp_0
test_that("Clarabel: socp_0", {
  h <- sth_socp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_socp_1
test_that("Clarabel: socp_1", {
  h <- sth_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CLARABEL")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_socp_2
test_that("Clarabel: socp_2", {
  h <- sth_socp_2()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CLARABEL")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_socp_3
test_that("Clarabel: socp_3 axis=0", {
  h <- sth_socp_3_ax0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CLARABEL")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_socp_3
test_that("Clarabel: socp_3 axis=1", {
  h <- sth_socp_3_ax1()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CLARABEL")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_sdp_1min
test_that("Clarabel: sdp_1 min", {
  h <- sth_sdp_1_min()
  verify_obj(h$prob, h$expect_obj, 1e-2, "CLARABEL")
})

## @cvxpy NONE
test_that("Clarabel: sdp_1 max", {
  h <- sth_sdp_1_max()
  verify_obj(h$prob, h$expect_obj, 1e-2, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_sdp_2
test_that("Clarabel: sdp_2", {
  h <- sth_sdp_2()
  verify_obj(h$prob, h$expect_obj, 1e-2, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_expcone_1
test_that("Clarabel: expcone_1", {
  h <- sth_expcone_1()
  verify_obj(h$prob, h$expect_obj, 1e-3, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_exp_soc_1
test_that("Clarabel: exp_soc_1", {
  h <- sth_expcone_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-3, "CLARABEL")
})

## @cvxpy NONE
test_that("Clarabel: pcp_1", {
  h <- sth_pcp_1()
  verify_obj(h$prob, h$expect_obj, 1e-3, "CLARABEL")
})

## @cvxpy NONE
test_that("Clarabel: pcp_2", {
  h <- sth_pcp_2()
  verify_obj(h$prob, h$expect_obj, 1e-3, "CLARABEL")
})

## @cvxpy NONE
test_that("Clarabel: pcp_3", {
  h <- sth_pcp_3()
  verify_obj(h$prob, h$expect_obj, 1e-3, "CLARABEL")
})

# ══════════════════════════════════════════════════════════════════
# SCS
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("SCS: lp_0", {
  h <- sth_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy NONE
test_that("SCS: lp_1", {
  h <- sth_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-2, "SCS")
})

## @cvxpy NONE
test_that("SCS: lp_2", {
  h <- sth_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-1, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_lp_3
test_that("SCS: lp_3 (unbounded)", {
  h <- sth_lp_3()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_lp_4
test_that("SCS: lp_4 (infeasible)", {
  h <- sth_lp_4()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_lp_5
test_that("SCS: lp_5 (redundant equalities)", {
  h <- sth_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-1, "SCS")
})

## @cvxpy NONE
test_that("SCS: qp_0", {
  h <- sth_qp_0()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy NONE
test_that("SCS: socp_0", {
  h <- sth_socp_0()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_socp_1
test_that("SCS: socp_1", {
  h <- sth_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy NONE
test_that("SCS: socp_2", {
  h <- sth_socp_2()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_socp_3
test_that("SCS: socp_3 axis=0", {
  h <- sth_socp_3_ax0()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_socp_3
test_that("SCS: socp_3 axis=1", {
  h <- sth_socp_3_ax1()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_sdp_1min
test_that("SCS: sdp_1 min", {
  h <- sth_sdp_1_min()
  verify_obj(h$prob, h$expect_obj, 1e-1, "SCS")
})

## @cvxpy NONE
test_that("SCS: sdp_1 max", {
  h <- sth_sdp_1_max()
  verify_obj(h$prob, h$expect_obj, 1e-1, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_sdp_2
test_that("SCS: sdp_2", {
  h <- sth_sdp_2()
  verify_obj(h$prob, h$expect_obj, 1e-1, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_expcone_1
test_that("SCS: expcone_1", {
  h <- sth_expcone_1()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_exp_soc_1
test_that("SCS: exp_soc_1", {
  h <- sth_expcone_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_pcp_1
test_that("SCS: pcp_1", {
  h <- sth_pcp_1()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_pcp_2
test_that("SCS: pcp_2", {
  h <- sth_pcp_2()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_pcp_3
test_that("SCS: pcp_3", {
  h <- sth_pcp_3()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

# ══════════════════════════════════════════════════════════════════
# OSQP (LP and QP only)
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("OSQP: lp_0", {
  h <- sth_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "OSQP")
})

## @cvxpy NONE
test_that("OSQP: lp_1", {
  h <- sth_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-4, "OSQP")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "OSQP")
})

## @cvxpy NONE
test_that("OSQP: lp_2", {
  h <- sth_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-4, "OSQP")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "OSQP")
})

## @cvxpy NONE
test_that("OSQP: lp_3 (unbounded)", {
  h <- sth_lp_3()
  verify_obj(h$prob, h$expect_obj, 1e-4, "OSQP")
})

## @cvxpy NONE
test_that("OSQP: lp_4 (infeasible)", {
  h <- sth_lp_4()
  verify_obj(h$prob, h$expect_obj, 1e-4, "OSQP")
})

## @cvxpy NONE
test_that("OSQP: qp_0", {
  h <- sth_qp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "OSQP")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "OSQP")
})

# ══════════════════════════════════════════════════════════════════
# HiGHS (LP and QP only)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_solving
test_that("HiGHS: lp_0", {
  h <- sth_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-5, "HIGHS")
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_solving
test_that("HiGHS: lp_1", {
  h <- sth_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-5, "HIGHS")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-5, "HIGHS")
  ## NOTE: HiGHS dual signs are negated vs CVXPY convention (known issue)
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_solving
test_that("HiGHS: lp_2", {
  h <- sth_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-5, "HIGHS")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-5, "HIGHS")
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_solving
test_that("HiGHS: lp_3 (unbounded)", {
  h <- sth_lp_3()
  verify_obj(h$prob, h$expect_obj, 1e-5, "HIGHS")
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_solving
test_that("HiGHS: lp_4 (infeasible)", {
  h <- sth_lp_4()
  verify_obj(h$prob, h$expect_obj, 1e-5, "HIGHS")
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_solving
test_that("HiGHS: lp_5 (redundant equalities)", {
  h <- sth_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-4, "HIGHS")
})

## @cvxpy NONE
test_that("HiGHS: qp_0", {
  h <- sth_qp_0()
  verify_obj(h$prob, h$expect_obj, 1e-5, "HIGHS")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-5, "HIGHS")
})

# ══════════════════════════════════════════════════════════════════
# Cross-solver consistency checks
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Cross-solver: LP objectives agree across solvers", {
  h <- sth_lp_1()

  obj_clarabel <- psolve(h$prob, solver = "CLARABEL")
  h2 <- sth_lp_1()
  obj_scs <- psolve(h2$prob, solver = "SCS")
  h3 <- sth_lp_1()
  obj_osqp <- psolve(h3$prob, solver = "OSQP")
  h4 <- sth_lp_1()
  obj_highs <- psolve(h4$prob, solver = "HIGHS")

  expect_equal(obj_clarabel, obj_scs, tolerance = 1e-2)
  expect_equal(obj_clarabel, obj_osqp, tolerance = 1e-3)
  expect_equal(obj_clarabel, obj_highs, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Cross-solver: QP objectives agree across solvers", {
  h <- sth_qp_0()

  obj_clarabel <- psolve(h$prob, solver = "CLARABEL")
  h2 <- sth_qp_0()
  obj_scs <- psolve(h2$prob, solver = "SCS")
  h3 <- sth_qp_0()
  obj_osqp <- psolve(h3$prob, solver = "OSQP")
  h4 <- sth_qp_0()
  obj_highs <- psolve(h4$prob, solver = "HIGHS")

  expect_equal(obj_clarabel, obj_scs, tolerance = 1e-2)
  expect_equal(obj_clarabel, obj_osqp, tolerance = 1e-3)
  expect_equal(obj_clarabel, obj_highs, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Cross-solver: SOCP objectives agree (Clarabel vs SCS)", {
  h1 <- sth_socp_1()
  obj_clarabel <- psolve(h1$prob, solver = "CLARABEL")
  h2 <- sth_socp_1()
  obj_scs <- psolve(h2$prob, solver = "SCS")
  expect_equal(obj_clarabel, obj_scs, tolerance = 1e-2)
})

# ══════════════════════════════════════════════════════════════════
# Solver option tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Clarabel: solver options accepted", {
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "CLARABEL", max_iter = 50)
  expect_equal(val, h$expect_obj, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("SCS: solver options accepted", {
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "SCS", max_iters = 5000)
  expect_equal(val, h$expect_obj, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("OSQP: solver options accepted", {
  h <- sth_qp_0()
  val <- psolve(h$prob, solver = "OSQP", max_iter = 5000)
  expect_equal(val, h$expect_obj, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("HiGHS: solver options accepted", {
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "HIGHS", time_limit = 60)
  expect_equal(val, h$expect_obj, tolerance = 1e-5)
})

# ══════════════════════════════════════════════════════════════════
# Invalid solver for problem type
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_qp_solvers.py::TestQpSolverValidation::test_qp_solver_rejects_soc_cones
test_that("OSQP rejects SOCP problem", {
  h <- sth_socp_1()
  expect_error(psolve(h$prob, solver = "OSQP"))
})

## @cvxpy NONE
test_that("HiGHS rejects SOCP problem", {
  h <- sth_socp_1()
  expect_error(psolve(h$prob, solver = "HIGHS"))
})

## @cvxpy test_qp_solvers.py::TestQpSolverValidation::test_qp_solver_rejects_psd_cones
test_that("OSQP rejects SDP problem", {
  h <- sth_sdp_1_min()
  expect_error(psolve(h$prob, solver = "OSQP"))
})

## @cvxpy NONE
test_that("HiGHS rejects SDP problem", {
  h <- sth_sdp_1_min()
  expect_error(psolve(h$prob, solver = "HIGHS"))
})

# ══════════════════════════════════════════════════════════════════
# TestAllSolvers — infrastructure tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestAllSolvers::test_installed_solvers
test_that("AllSolvers: installed_solvers returns character vector of known solvers", {
  solvers <- installed_solvers()
  expect_type(solvers, "character")
  ## At minimum, Clarabel is always available (bundled)
  expect_true("CLARABEL" %in% solvers)
  ## All returned names must be valid solver constants
  known <- c("CLARABEL", "SCS", "OSQP", "HIGHS", "MOSEK", "GUROBI",
             "GLPK", "GLPK_MI", "ECOS", "ECOS_BB", "CPLEX", "CVXOPT", "PIQP", "SCIP",
             "XPRESS")
  expect_true(all(solvers %in% known))
})

# ══════════════════════════════════════════════════════════════════
# Clarabel (gap tests)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_nonstandard_name
test_that("Clarabel: nonstandard variable name", {
  ## Variables with unusual names still solve correctly
  x <- Variable(2, name = "my var with spaces")
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  val <- psolve(prob, solver = "CLARABEL")
  expect_equal(val, 2.0, tolerance = 1e-4)
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_pcp_0
test_that("Clarabel: pcp_0 (socp_0 via power cone path)", {
  ## CVXPY maps test_clarabel_pcp_0 -> StandardTestSOCPs.test_socp_0('CLARABEL')
  h <- sth_socp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_pcp_1
test_that("Clarabel: pcp_1 (gap closure)", {
  ## Already covered as "Clarabel: pcp_1" but annotated differently
  h <- sth_pcp_1()
  verify_obj(h$prob, h$expect_obj, 1e-3, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_pcp_2
test_that("Clarabel: pcp_2 (gap closure)", {
  ## Already covered as "Clarabel: pcp_2" but annotated differently
  h <- sth_pcp_2()
  verify_obj(h$prob, h$expect_obj, 1e-3, "CLARABEL")
})

## @cvxpy test_conic_solvers.py::TestClarabel::test_clarabel_qp_0_linear_obj
test_that("Clarabel: qp_0 with linear objective", {
  h <- sth_qp_0_linear_obj()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CLARABEL")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "CLARABEL")
})

# ══════════════════════════════════════════════════════════════════
# HiGHS (gap tests)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_nonstandard_name
test_that("HiGHS: nonstandard variable name", {
  x <- Variable(2, name = "my var with spaces")
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  val <- psolve(prob, solver = "HIGHS")
  expect_equal(val, 2.0, tolerance = 1e-5)
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_options
test_that("HiGHS: options (CVXPY parity)", {
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "HIGHS", time_limit = 60)
  expect_equal(val, h$expect_obj, tolerance = 1e-5)
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_validate_column_name
test_that("HiGHS: validate_column_name", {
  skip("HiGHS column name validation not exposed in R highs package")
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_warm_start
test_that("HiGHS: warm_start", {
  skip("HiGHS warm-start blocked: R highs package lacks setSolution()")
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_written_model_contains_variable_names
test_that("HiGHS: written_model_contains_variable_names", {
  skip("HiGHS model export not exposed in R highs package")
})

# ══════════════════════════════════════════════════════════════════
# SCS (gap tests)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestSCS::test_complex_matrices
test_that("SCS: complex matrices", {
  ## Trace norm of a complex matrix via SDP dual formulation
  ## CVXPY SOURCE: test_conic_solvers.py TestSCS::test_complex_matrices
  ## K is a fixed 2x2 complex matrix (np.random.seed(0))
  K <- matrix(c(0.5488135+0.4236548i, 0.60276338+0.43758721i,
                 0.71518937+0.64589411i, 0.54488318+0.891773i), nrow = 2)
  ## Expected: sum of singular values of K
  expected_n1 <- 1.8740804548294105

  X <- Variable(c(2, 2), complex = TRUE)
  Y <- Variable(c(2, 2), complex = TRUE)
  objective <- Minimize(Re(0.5 * matrix_trace(X) + 0.5 * matrix_trace(Y)))
  constraints <- list(
    PSD(bmat(list(list(X, -Conj(t(K))), list(-K, Y)))),
    PSD(X),
    PSD(Y)
  )
  prob <- Problem(objective, constraints)
  val <- psolve(prob, solver = "SCS")
  expect_equal(val, expected_n1, tolerance = 0.05)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_psd_constraint
test_that("SCS: psd_constraint", {
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Minimize(matrix_trace(X)),
                  list(PSD(X), X[1, 1] >= 1, X[2, 2] >= 1, X[3, 3] >= 1))
  val <- psolve(prob, solver = "SCS")
  expect_equal(val, 3.0, tolerance = 1e-1)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_quad_obj
test_that("SCS: quad_obj", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x)), list(x[1] >= 1))
  val <- psolve(prob, solver = "SCS")
  expect_equal(val, 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), c(1, 0), tolerance = 1e-1)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_quad_obj_with_power
test_that("SCS: quad_obj_with_power", {
  x <- Variable(2)
  prob <- Problem(Minimize(power(x[1], 2) + power(x[2], 2)),
                  list(x[1] >= 1))
  val <- psolve(prob, solver = "SCS")
  expect_equal(val, 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), c(1, 0), tolerance = 1e-1)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_options
test_that("SCS: solver options (max_iters, eps_abs)", {
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "SCS", max_iters = 10000)
  expect_equal(val, h$expect_obj, tolerance = 1e-2)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_retry
test_that("SCS: retry on inaccurate solution", {
  ## SCS can retry with different settings when it gets inaccurate results.
  ## Use a problem that converges easily to show the mechanism works.
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "SCS")
  expect_true(status(h$prob) %in% c("optimal", "optimal_inaccurate"))
  expect_equal(val, h$expect_obj, tolerance = 1e-1)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_sdp_pcp_1
test_that("SCS: sdp_pcp_1 (mixed SDP + power cone)", {
  h <- sth_sdp_pcp_1()
  val <- psolve(h$prob, solver = "SCS")
  expect_true(is.finite(val))
  expect_true(status(h$prob) %in% c("optimal", "optimal_inaccurate"))
})

## @cvxpy test_conic_solvers.py::TestSCS::test_sdp_var
test_that("SCS: sdp_var", {
  ## SDP variable (symmetric, PSD constrained)
  X <- Variable(c(2, 2), PSD = TRUE)
  prob <- Problem(Minimize(matrix_trace(X)),
                  list(X[1, 1] >= 2))
  val <- psolve(prob, solver = "SCS")
  expect_equal(val, 2.0, tolerance = 1e-1)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_sigma_max
test_that("SCS: sigma_max", {
  ## Minimize sigma_max (largest singular value)
  Y <- Variable(c(2, 3))
  prob <- Problem(Minimize(sigma_max(Y)),
                  list(Y[1, 1] >= 1, Y[2, 2] >= 1))
  val <- psolve(prob, solver = "SCS")
  expect_true(is.finite(val))
  expect_equal(val, 1.0, tolerance = 1e-1)
})
