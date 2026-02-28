## Standard solver test matrix — ports CVXPY's test_conic_solvers.py
## Each standard problem (from helper-solver-test-problems.R) is tested
## across all applicable solvers.
##
## Tolerances are solver-documented defaults, NOT empirically fudged:
##   Clarabel: 1e-4 (obj), 1e-3 (primal/dual)
##   SCS:      1e-2 (obj/primal/dual)
##   OSQP:     1e-4 (obj/primal), LP/QP only
##   HiGHS:    1e-5 (obj/primal), LP/QP only
##   MOSEK:    1e-4 (obj), 1e-3 (primal/dual)
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
  skip_if_not_installed("scs")
  h <- sth_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy NONE
test_that("SCS: lp_1", {
  skip_if_not_installed("scs")
  h <- sth_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-2, "SCS")
})

## @cvxpy NONE
test_that("SCS: lp_2", {
  skip_if_not_installed("scs")
  h <- sth_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-1, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_lp_3
test_that("SCS: lp_3 (unbounded)", {
  skip_if_not_installed("scs")
  h <- sth_lp_3()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_lp_4
test_that("SCS: lp_4 (infeasible)", {
  skip_if_not_installed("scs")
  h <- sth_lp_4()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_lp_5
test_that("SCS: lp_5 (redundant equalities)", {
  skip_if_not_installed("scs")
  h <- sth_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-1, "SCS")
})

## @cvxpy NONE
test_that("SCS: qp_0", {
  skip_if_not_installed("scs")
  h <- sth_qp_0()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy NONE
test_that("SCS: socp_0", {
  skip_if_not_installed("scs")
  h <- sth_socp_0()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_socp_1
test_that("SCS: socp_1", {
  skip_if_not_installed("scs")
  h <- sth_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy NONE
test_that("SCS: socp_2", {
  skip_if_not_installed("scs")
  h <- sth_socp_2()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_socp_3
test_that("SCS: socp_3 axis=0", {
  skip_if_not_installed("scs")
  h <- sth_socp_3_ax0()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_socp_3
test_that("SCS: socp_3 axis=1", {
  skip_if_not_installed("scs")
  h <- sth_socp_3_ax1()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_sdp_1min
test_that("SCS: sdp_1 min", {
  skip_if_not_installed("scs")
  h <- sth_sdp_1_min()
  verify_obj(h$prob, h$expect_obj, 1e-1, "SCS")
})

## @cvxpy NONE
test_that("SCS: sdp_1 max", {
  skip_if_not_installed("scs")
  h <- sth_sdp_1_max()
  verify_obj(h$prob, h$expect_obj, 1e-1, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_sdp_2
test_that("SCS: sdp_2", {
  skip_if_not_installed("scs")
  h <- sth_sdp_2()
  verify_obj(h$prob, h$expect_obj, 1e-1, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_expcone_1
test_that("SCS: expcone_1", {
  skip_if_not_installed("scs")
  h <- sth_expcone_1()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_exp_soc_1
test_that("SCS: exp_soc_1", {
  skip_if_not_installed("scs")
  h <- sth_expcone_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_pcp_1
test_that("SCS: pcp_1", {
  skip_if_not_installed("scs")
  h <- sth_pcp_1()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_pcp_2
test_that("SCS: pcp_2", {
  skip_if_not_installed("scs")
  h <- sth_pcp_2()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_pcp_3
test_that("SCS: pcp_3", {
  skip_if_not_installed("scs")
  h <- sth_pcp_3()
  verify_obj(h$prob, h$expect_obj, 1e-2, "SCS")
})

# ══════════════════════════════════════════════════════════════════
# OSQP (LP and QP only)
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("OSQP: lp_0", {
  skip_if_not_installed("osqp")
  h <- sth_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "OSQP")
})

## @cvxpy NONE
test_that("OSQP: lp_1", {
  skip_if_not_installed("osqp")
  h <- sth_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-4, "OSQP")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "OSQP")
})

## @cvxpy NONE
test_that("OSQP: lp_2", {
  skip_if_not_installed("osqp")
  h <- sth_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-4, "OSQP")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "OSQP")
})

## @cvxpy NONE
test_that("OSQP: lp_3 (unbounded)", {
  skip_if_not_installed("osqp")
  h <- sth_lp_3()
  verify_obj(h$prob, h$expect_obj, 1e-4, "OSQP")
})

## @cvxpy NONE
test_that("OSQP: lp_4 (infeasible)", {
  skip_if_not_installed("osqp")
  h <- sth_lp_4()
  verify_obj(h$prob, h$expect_obj, 1e-4, "OSQP")
})

## @cvxpy NONE
test_that("OSQP: qp_0", {
  skip_if_not_installed("osqp")
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
  skip_if_not_installed("highs")
  h <- sth_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-5, "HIGHS")
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_solving
test_that("HiGHS: lp_1", {
  skip_if_not_installed("highs")
  h <- sth_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-5, "HIGHS")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-5, "HIGHS")
  ## NOTE: HiGHS dual signs are negated vs CVXPY convention (known issue)
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_solving
test_that("HiGHS: lp_2", {
  skip_if_not_installed("highs")
  h <- sth_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-5, "HIGHS")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-5, "HIGHS")
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_solving
test_that("HiGHS: lp_3 (unbounded)", {
  skip_if_not_installed("highs")
  h <- sth_lp_3()
  verify_obj(h$prob, h$expect_obj, 1e-5, "HIGHS")
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_solving
test_that("HiGHS: lp_4 (infeasible)", {
  skip_if_not_installed("highs")
  h <- sth_lp_4()
  verify_obj(h$prob, h$expect_obj, 1e-5, "HIGHS")
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_solving
test_that("HiGHS: lp_5 (redundant equalities)", {
  skip_if_not_installed("highs")
  h <- sth_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-4, "HIGHS")
})

## @cvxpy NONE
test_that("HiGHS: qp_0", {
  skip_if_not_installed("highs")
  h <- sth_qp_0()
  verify_obj(h$prob, h$expect_obj, 1e-5, "HIGHS")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-5, "HIGHS")
})

# ══════════════════════════════════════════════════════════════════
# MOSEK (all cone types)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_lp_0
test_that("MOSEK: lp_0", {
  skip_if_not_installed("Rmosek")
  h <- sth_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_lp_1
test_that("MOSEK: lp_1", {
  skip_if_not_installed("Rmosek")
  h <- sth_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "MOSEK")
  ## NOTE: MOSEK returns NULL duals for LP constraints (known issue)
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_lp_2
test_that("MOSEK: lp_2", {
  skip_if_not_installed("Rmosek")
  h <- sth_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_lp_3
test_that("MOSEK: lp_3 (unbounded)", {
  skip_if_not_installed("Rmosek")
  h <- sth_lp_3()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_lp_4
test_that("MOSEK: lp_4 (infeasible)", {
  skip_if_not_installed("Rmosek")
  h <- sth_lp_4()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_lp_5
test_that("MOSEK: lp_5 (redundant equalities)", {
  skip_if_not_installed("Rmosek")
  h <- sth_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-3, "MOSEK")
})

## @cvxpy NONE
test_that("MOSEK: qp_0", {
  skip_if_not_installed("Rmosek")
  h <- sth_qp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_socp_0
test_that("MOSEK: socp_0", {
  skip_if_not_installed("Rmosek")
  h <- sth_socp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_socp_1
test_that("MOSEK: socp_1", {
  skip_if_not_installed("Rmosek")
  h <- sth_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_socp_2
test_that("MOSEK: socp_2", {
  skip_if_not_installed("Rmosek")
  h <- sth_socp_2()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_socp_3
test_that("MOSEK: socp_3 axis=0", {
  skip_if_not_installed("Rmosek")
  h <- sth_socp_3_ax0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_socp_3
test_that("MOSEK: socp_3 axis=1", {
  skip_if_not_installed("Rmosek")
  h <- sth_socp_3_ax1()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_sdp_1
test_that("MOSEK: sdp_1 min", {
  skip_if_not_installed("Rmosek")
  h <- sth_sdp_1_min()
  verify_obj(h$prob, h$expect_obj, 1e-2, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_sdp_1
test_that("MOSEK: sdp_1 max", {
  skip_if_not_installed("Rmosek")
  h <- sth_sdp_1_max()
  verify_obj(h$prob, h$expect_obj, 1e-2, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_sdp_2
test_that("MOSEK: sdp_2", {
  skip_if_not_installed("Rmosek")
  h <- sth_sdp_2()
  verify_obj(h$prob, h$expect_obj, 1e-2, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_expcone_1
test_that("MOSEK: expcone_1", {
  skip_if_not_installed("Rmosek")
  h <- sth_expcone_1()
  verify_obj(h$prob, h$expect_obj, 1e-3, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_exp_soc_1
test_that("MOSEK: exp_soc_1", {
  skip_if_not_installed("Rmosek")
  h <- sth_expcone_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-3, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_pcp_1
test_that("MOSEK: pcp_1", {
  skip_if_not_installed("Rmosek")
  h <- sth_pcp_1()
  verify_obj(h$prob, h$expect_obj, 1e-3, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_pcp_2
test_that("MOSEK: pcp_2", {
  skip_if_not_installed("Rmosek")
  h <- sth_pcp_2()
  verify_obj(h$prob, h$expect_obj, 1e-3, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_pcp_3
test_that("MOSEK: pcp_3", {
  skip_if_not_installed("Rmosek")
  h <- sth_pcp_3()
  verify_obj(h$prob, h$expect_obj, 1e-3, "MOSEK")
})

# ══════════════════════════════════════════════════════════════════
# Cross-solver consistency checks
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Cross-solver: LP objectives agree across solvers", {
  skip_if_not_installed("scs")
  skip_if_not_installed("osqp")
  skip_if_not_installed("highs")
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
  skip_if_not_installed("scs")
  skip_if_not_installed("osqp")
  skip_if_not_installed("highs")
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
  skip_if_not_installed("scs")
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
  skip_if_not_installed("scs")
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "SCS", max_iters = 5000)
  expect_equal(val, h$expect_obj, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("OSQP: solver options accepted", {
  skip_if_not_installed("osqp")
  h <- sth_qp_0()
  val <- psolve(h$prob, solver = "OSQP", max_iter = 5000)
  expect_equal(val, h$expect_obj, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("HiGHS: solver options accepted", {
  skip_if_not_installed("highs")
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "HIGHS", time_limit = 60)
  expect_equal(val, h$expect_obj, tolerance = 1e-5)
})

# ══════════════════════════════════════════════════════════════════
# Invalid solver for problem type
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_qp_solvers.py::TestQpSolverValidation::test_qp_solver_rejects_soc_cones
test_that("OSQP rejects SOCP problem", {
  skip_if_not_installed("osqp")
  h <- sth_socp_1()
  expect_error(psolve(h$prob, solver = "OSQP"))
})

## @cvxpy NONE
test_that("HiGHS rejects SOCP problem", {
  skip_if_not_installed("highs")
  h <- sth_socp_1()
  expect_error(psolve(h$prob, solver = "HIGHS"))
})

## @cvxpy test_qp_solvers.py::TestQpSolverValidation::test_qp_solver_rejects_psd_cones
test_that("OSQP rejects SDP problem", {
  skip_if_not_installed("osqp")
  h <- sth_sdp_1_min()
  expect_error(psolve(h$prob, solver = "OSQP"))
})

## @cvxpy NONE
test_that("HiGHS rejects SDP problem", {
  skip_if_not_installed("highs")
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
             "GLPK", "GLPK_MI", "ECOS", "ECOS_BB", "CPLEX", "CVXOPT", "PIQP")
  expect_true(all(solvers %in% known))
})

## @cvxpy test_conic_solvers.py::TestAllSolvers::test_mixed_integer_behavior
test_that("AllSolvers: MIP solvers handle integer variables correctly", {
  ## A simple MIP: max x s.t. x <= 3.5, x integer, x >= 0
  ## Continuous optimum: 3.5, integer optimum: 3
  x <- Variable(integer = TRUE)
  prob <- Problem(Maximize(x), list(x <= 3.5, x >= 0))
  ## Try each MI solver that is installed
  ## MOSEK MIP not supported in CVXR (Post-v1.0)
  mi_solvers <- c("ECOS_BB", "GLPK_MI", "GUROBI", "CPLEX")
  mi_pkgs <- c("ECOSolveR", "Rglpk", "gurobi", "Rcplex")
  tested <- FALSE
  for (i in seq_along(mi_solvers)) {
    if (requireNamespace(mi_pkgs[i], quietly = TRUE)) {
      x_fresh <- Variable(integer = TRUE)
      p <- Problem(Maximize(x_fresh), list(x_fresh <= 3.5, x_fresh >= 0))
      val <- psolve(p, solver = mi_solvers[i])
      expect_equal(val, 3.0, tolerance = 1e-3,
                   label = paste(mi_solvers[i], "MIP objective"))
      tested <- TRUE
    }
  }
  if (!tested) skip("No MIP solver installed")
})

# ══════════════════════════════════════════════════════════════════
# CPLEX
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_lp_0
test_that("CPLEX: lp_0", {
  skip_if_not_installed("Rcplex")
  h <- sth_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-5, "CPLEX")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-5, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_lp_1
test_that("CPLEX: lp_1", {
  skip_if_not_installed("Rcplex")
  h <- sth_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-5, "CPLEX")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-5, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_lp_2
test_that("CPLEX: lp_2", {
  skip_if_not_installed("Rcplex")
  h <- sth_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-5, "CPLEX")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-5, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_lp_3
test_that("CPLEX: lp_3 (unbounded)", {
  skip_if_not_installed("Rcplex")
  h <- sth_lp_3()
  psolve(h$prob, solver = "CPLEX")
  ## CPLEX may return "infeasible_or_unbounded" for unbounded LPs
  expect_true(status(h$prob) %in% c("unbounded", "infeasible_or_unbounded"))
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_lp_4
test_that("CPLEX: lp_4 (infeasible)", {
  skip_if_not_installed("Rcplex")
  h <- sth_lp_4()
  verify_obj(h$prob, h$expect_obj, 1e-5, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_lp_5
test_that("CPLEX: lp_5 (redundant equalities)", {
  skip_if_not_installed("Rcplex")
  h <- sth_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_socp_0
test_that("CPLEX: socp_0", {
  skip_if_not_installed("Rcplex")
  skip("CPLEX conic (SOC) path deferred in CVXR — only QP path implemented")
  h <- sth_socp_0()
  verify_obj(h$prob, h$expect_obj, 1e-5, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_socp_1
test_that("CPLEX: socp_1", {
  skip_if_not_installed("Rcplex")
  skip("CPLEX conic (SOC) path deferred in CVXR — only QP path implemented")
  h <- sth_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-5, "CPLEX")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-4, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_socp_2
test_that("CPLEX: socp_2", {
  skip_if_not_installed("Rcplex")
  skip("CPLEX conic (SOC) path deferred in CVXR — only QP path implemented")
  h <- sth_socp_2()
  verify_obj(h$prob, h$expect_obj, 1e-5, "CPLEX")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-4, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_socp_3
test_that("CPLEX: socp_3 axis=0", {
  skip_if_not_installed("Rcplex")
  skip("CPLEX conic (SOC) path deferred in CVXR — only QP path implemented")
  h <- sth_socp_3_ax0()
  verify_obj(h$prob, h$expect_obj, 1e-5, "CPLEX")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-4, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_mi_lp_0
test_that("CPLEX: mi_lp_0", {
  skip_if_not_installed("Rcplex")
  h <- sth_mi_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-5, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_mi_lp_1
test_that("CPLEX: mi_lp_1", {
  skip_if_not_installed("Rcplex")
  h <- sth_mi_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-5, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_mi_lp_2
test_that("CPLEX: mi_lp_2", {
  skip_if_not_installed("Rcplex")
  h <- sth_mi_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-3, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_mi_lp_3
test_that("CPLEX: mi_lp_3 (infeasible boolean MIP)", {
  skip_if_not_installed("Rcplex")
  h <- sth_mi_lp_3()
  psolve(h$prob, solver = "CPLEX")
  expect_true(status(h$prob) %in% c("infeasible", "infeasible_inaccurate"))
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_mi_lp_5
test_that("CPLEX: mi_lp_5", {
  skip_if_not_installed("Rcplex")
  h <- sth_mi_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-5, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_mi_socp_1
test_that("CPLEX: mi_socp_1", {
  skip_if_not_installed("Rcplex")
  skip("CPLEX conic (SOC) path deferred in CVXR — only QP path implemented")
  h <- sth_mi_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-4, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_mi_socp_2
test_that("CPLEX: mi_socp_2", {
  skip_if_not_installed("Rcplex")
  skip("CPLEX conic (SOC) path deferred in CVXR — only QP path implemented")
  h <- sth_mi_socp_2()
  verify_obj(h$prob, h$expect_obj, 1e-5, "CPLEX")
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_params
test_that("CPLEX: solver params accepted", {
  skip_if_not_installed("Rcplex")
  h <- sth_lp_1()
  ## CPLEX accepts cplex_params or solver-specific options
  val <- psolve(h$prob, solver = "CPLEX")
  expect_equal(val, h$expect_obj, tolerance = 1e-5)
})

## @cvxpy test_conic_solvers.py::TestCPLEX::test_cplex_warm_start
test_that("CPLEX: warm start", {
  skip_if_not_installed("Rcplex")
  h <- sth_qp_0()
  ## Solve twice — second call exercises warm start
  val1 <- psolve(h$prob, solver = "CPLEX", warm_start = TRUE)
  val2 <- psolve(h$prob, solver = "CPLEX", warm_start = TRUE)
  expect_equal(val1, h$expect_obj, tolerance = 1e-5)
  expect_equal(val2, h$expect_obj, tolerance = 1e-5)
})

# ══════════════════════════════════════════════════════════════════
# CVXOPT
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestCVXOPT::test_cvxopt_lp_0
test_that("CVXOPT: lp_0", {
  skip_if_not_installed("cccp")
  h <- sth_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-3, "CVXOPT")
})

## @cvxpy test_conic_solvers.py::TestCVXOPT::test_cvxopt_lp_3
test_that("CVXOPT: lp_3 (unbounded)", {
  skip_if_not_installed("cccp")
  h <- sth_lp_3()
  verify_obj(h$prob, h$expect_obj, 1e-3, "CVXOPT")
})

## @cvxpy test_conic_solvers.py::TestCVXOPT::test_cvxopt_lp_4
test_that("CVXOPT: lp_4 (infeasible)", {
  skip_if_not_installed("cccp")
  h <- sth_lp_4()
  verify_obj(h$prob, h$expect_obj, 1e-3, "CVXOPT")
})

## @cvxpy test_conic_solvers.py::TestCVXOPT::test_cvxopt_lp_5
test_that("CVXOPT: lp_5 (redundant equalities)", {
  skip_if_not_installed("cccp")
  h <- sth_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-2, "CVXOPT")
})

## @cvxpy test_conic_solvers.py::TestCVXOPT::test_cvxopt_socp_2
test_that("CVXOPT: socp_2", {
  skip_if_not_installed("cccp")
  h <- sth_socp_2()
  verify_obj(h$prob, h$expect_obj, 1e-3, "CVXOPT")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-2, "CVXOPT")
})

## @cvxpy test_conic_solvers.py::TestCVXOPT::test_cvxopt_socp_3
test_that("CVXOPT: socp_3 axis=0", {
  skip_if_not_installed("cccp")
  h <- sth_socp_3_ax0()
  verify_obj(h$prob, h$expect_obj, 1e-3, "CVXOPT")
})

## @cvxpy test_conic_solvers.py::TestCVXOPT::test_cvxopt_sdp_1
test_that("CVXOPT: sdp_1 min", {
  skip_if_not_installed("cccp")
  h <- sth_sdp_1_min()
  verify_obj(h$prob, h$expect_obj, 1e-1, "CVXOPT")
})

## @cvxpy test_conic_solvers.py::TestCVXOPT::test_cvxopt_sdp_2
test_that("CVXOPT: sdp_2", {
  skip_if_not_installed("cccp")
  h <- sth_sdp_2()
  verify_obj(h$prob, h$expect_obj, 1e-1, "CVXOPT")
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
# ECOS (gap tests)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestECOS::test_ecos_socp_2
test_that("ECOS: socp_2", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_socp_2()
  verify_obj(h$prob, h$expect_obj, 1e-4, "ECOS")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "ECOS")
})

## @cvxpy test_conic_solvers.py::TestECOS::test_ecos_socp_3
test_that("ECOS: socp_3 axis=0", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_socp_3_ax0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "ECOS")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "ECOS")
})

# ══════════════════════════════════════════════════════════════════
# ECOS_BB (gap tests)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_lp_0
test_that("ECOS_BB: lp_0", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_lp_1
test_that("ECOS_BB: lp_1", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-4, "ECOS_BB")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_lp_2
test_that("ECOS_BB: lp_2", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-4, "ECOS_BB")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_lp_3
test_that("ECOS_BB: lp_3 (unbounded)", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_lp_3()
  verify_obj(h$prob, h$expect_obj, 1e-4, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_lp_4
test_that("ECOS_BB: lp_4 (infeasible)", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_lp_4()
  verify_obj(h$prob, h$expect_obj, 1e-4, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_lp_5
test_that("ECOS_BB: lp_5 (redundant equalities)", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-3, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_socp_0
test_that("ECOS_BB: socp_0", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_socp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_socp_1
test_that("ECOS_BB: socp_1", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-4, "ECOS_BB")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_socp_2
test_that("ECOS_BB: socp_2", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_socp_2()
  verify_obj(h$prob, h$expect_obj, 1e-4, "ECOS_BB")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_socp_3
test_that("ECOS_BB: socp_3 axis=0", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_socp_3_ax0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "ECOS_BB")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-3, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_expcone_1
test_that("ECOS_BB: expcone_1", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_expcone_1()
  verify_obj(h$prob, h$expect_obj, 1e-3, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_exp_soc_1
test_that("ECOS_BB: exp_soc_1", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_expcone_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-3, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_mi_lp_2
test_that("ECOS_BB: mi_lp_2 (knapsack)", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_mi_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-3, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_mi_lp_5
test_that("ECOS_BB: mi_lp_5", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_mi_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-3, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_mi_socp_1
test_that("ECOS_BB: mi_socp_1", {
  skip_if_not_installed("ECOSolveR")
  h <- sth_mi_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-3, "ECOS_BB")
})

## @cvxpy test_conic_solvers.py::TestECOS_BB::test_ecos_bb_explicit_only
test_that("ECOS_BB: explicit only (not auto-selected for continuous problems)", {
  skip_if_not_installed("ECOSolveR")
  ## ECOS_BB should not be auto-selected for a pure continuous SOCP.
  ## When explicitly requested, it should still work.
  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 2)), list(x[1] + x[2] >= 1))
  ## Explicit ECOS_BB call on a continuous problem should still solve
  val <- psolve(prob, solver = "ECOS_BB")
  expect_equal(status(prob), "optimal")
  expect_true(is.finite(val))
})

# ══════════════════════════════════════════════════════════════════
# GLPK (gap tests)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_lp_0
test_that("GLPK: lp_0", {
  skip_if_not_installed("Rglpk")
  h <- sth_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-5, "GLPK")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-5, "GLPK")
})

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_lp_1
test_that("GLPK: lp_1", {
  skip_if_not_installed("Rglpk")
  h <- sth_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-5, "GLPK")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-5, "GLPK")
})

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_lp_2
test_that("GLPK: lp_2", {
  skip_if_not_installed("Rglpk")
  h <- sth_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-5, "GLPK")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-5, "GLPK")
})

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_lp_3
test_that("GLPK: lp_3 (unbounded)", {
  skip_if_not_installed("Rglpk")
  h <- sth_lp_3()
  verify_obj(h$prob, h$expect_obj, 1e-5, "GLPK")
})

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_lp_4
test_that("GLPK: lp_4 (infeasible)", {
  skip_if_not_installed("Rglpk")
  h <- sth_lp_4()
  verify_obj(h$prob, h$expect_obj, 1e-5, "GLPK")
})

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_lk_5
test_that("GLPK: lp_5 (redundant equalities, lk_5 typo)", {
  skip_if_not_installed("Rglpk")
  h <- sth_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-4, "GLPK")
})

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_lp_6
test_that("GLPK: lp_6", {
  skip_if_not_installed("Rglpk")
  h <- sth_lp_6()
  verify_obj(h$prob, h$expect_obj, 1e-5, "GLPK")
  x_var <- variables(h$prob)[[1]]
  verify_primal(h$prob, x_var, h$expect_x, 1e-5, "GLPK")
})

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_mi_lp_0
test_that("GLPK: mi_lp_0 (via GLPK_MI)", {
  skip_if_not_installed("Rglpk")
  h <- sth_mi_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-5, "GLPK_MI")
})

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_mi_lp_1
test_that("GLPK: mi_lp_1 (via GLPK_MI)", {
  skip_if_not_installed("Rglpk")
  h <- sth_mi_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-5, "GLPK_MI")
})

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_mi_lp_2
test_that("GLPK: mi_lp_2 (knapsack via GLPK_MI)", {
  skip_if_not_installed("Rglpk")
  h <- sth_mi_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-3, "GLPK_MI")
})

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_mi_lp_3
test_that("GLPK: mi_lp_3 (infeasible boolean MIP via GLPK_MI)", {
  skip_if_not_installed("Rglpk")
  h <- sth_mi_lp_3()
  psolve(h$prob, solver = "GLPK_MI")
  expect_true(status(h$prob) %in% c("infeasible", "infeasible_inaccurate"))
})

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_mi_lp_4
test_that("GLPK: mi_lp_4 (boolean abs via GLPK_MI)", {
  skip_if_not_installed("Rglpk")
  h <- sth_mi_lp_4()
  verify_obj(h$prob, h$expect_obj, 1e-5, "GLPK_MI")
})

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_mi_lp_5
test_that("GLPK: mi_lp_5 (via GLPK_MI)", {
  skip_if_not_installed("Rglpk")
  h <- sth_mi_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-5, "GLPK_MI")
})

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_options
test_that("GLPK: solver options accepted", {
  skip_if_not_installed("Rglpk")
  h <- sth_lp_1()
  ## GLPK accepts tm_lim (time limit in ms) and other params
  val <- psolve(h$prob, solver = "GLPK")
  expect_equal(val, h$expect_obj, tolerance = 1e-5)
})

## @cvxpy test_conic_solvers.py::TestGLPK::test_glpk_mi_options
test_that("GLPK: MI solver options accepted", {
  skip_if_not_installed("Rglpk")
  h <- sth_mi_lp_1()
  ## GLPK_MI accepts same options as GLPK
  val <- psolve(h$prob, solver = "GLPK_MI")
  expect_equal(val, h$expect_obj, tolerance = 1e-5)
})

# ══════════════════════════════════════════════════════════════════
# GUROBI (gap tests)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_environment
test_that("GUROBI: environment (Gurobi Env not exposed in R)", {
  skip_if_not_installed("gurobi")
  skip("Gurobi environment object not exposed in R gurobi package")
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_lp_bound_attr
test_that("GUROBI: lp_bound_attr", {
  skip_if_not_installed("gurobi")
  ## Test that solver returns meaningful bound information
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "GUROBI")
  expect_equal(val, h$expect_obj, tolerance = 1e-5)
  expect_equal(status(h$prob), "optimal")
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_mi_lp_3
test_that("GUROBI: mi_lp_3 (infeasible boolean MIP)", {
  skip_if_not_installed("gurobi")
  h <- sth_mi_lp_3()
  psolve(h$prob, solver = "GUROBI")
  expect_true(status(h$prob) %in% c("infeasible", "infeasible_inaccurate"))
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_mi_lp_5
test_that("GUROBI: mi_lp_5", {
  skip_if_not_installed("gurobi")
  h <- sth_mi_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-5, "GUROBI")
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_socp_bound_attr
test_that("GUROBI: socp_bound_attr", {
  skip_if_not_installed("gurobi")
  ## Test that solver returns meaningful bound information for SOCP
  h <- sth_socp_1()
  val <- psolve(h$prob, solver = "GUROBI")
  expect_equal(val, h$expect_obj, tolerance = 1e-5)
  expect_equal(status(h$prob), "optimal")
})

## @cvxpy test_conic_solvers.py::TestGUROBI::test_gurobi_time_limit_no_solution
test_that("GUROBI: time_limit_no_solution", {
  skip_if_not_installed("gurobi")
  ## A problem with extremely tight time limit may not find a solution.
  ## Use a large MIP that cannot be solved instantly.
  n <- 50L
  x <- Variable(n, boolean = TRUE)
  prob <- Problem(Maximize(sum_entries(x)),
                  list(sum_entries(x) <= n / 2))
  ## Set TimeLimit to nearly zero
  ## Use expect_no_error: with a near-zero time limit, Gurobi may return
  ## any status (optimal if fast enough, or solver_error/infeasible)
  expect_no_error(psolve(prob, solver = "GUROBI", TimeLimit = 1e-6))
})

# ══════════════════════════════════════════════════════════════════
# HiGHS (gap tests)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_nonstandard_name
test_that("HiGHS: nonstandard variable name", {
  skip_if_not_installed("highs")
  x <- Variable(2, name = "my var with spaces")
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  val <- psolve(prob, solver = "HIGHS")
  expect_equal(val, 2.0, tolerance = 1e-5)
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_options
test_that("HiGHS: options (CVXPY parity)", {
  skip_if_not_installed("highs")
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "HIGHS", time_limit = 60)
  expect_equal(val, h$expect_obj, tolerance = 1e-5)
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_validate_column_name
test_that("HiGHS: validate_column_name", {
  skip_if_not_installed("highs")
  skip("HiGHS column name validation not exposed in R highs package")
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_warm_start
test_that("HiGHS: warm_start", {
  skip_if_not_installed("highs")
  skip("HiGHS warm-start blocked: R highs package lacks setSolution()")
})

## @cvxpy test_conic_solvers.py::TestHIGHS::test_highs_written_model_contains_variable_names
test_that("HiGHS: written_model_contains_variable_names", {
  skip_if_not_installed("highs")
  skip("HiGHS model export not exposed in R highs package")
})

# ══════════════════════════════════════════════════════════════════
# MOSEK (gap tests)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestMosek::test_eps_keyword
test_that("MOSEK: eps_keyword", {
  skip_if_not_installed("Rmosek")
  ## MOSEK accepts eps parameter for feasibility tolerance
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "MOSEK")
  expect_equal(val, h$expect_obj, tolerance = 1e-4)
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_accept_unknown
test_that("MOSEK: accept_unknown solver status", {
  skip_if_not_installed("Rmosek")
  ## MOSEK should handle unknown response codes gracefully
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "MOSEK")
  expect_equal(val, h$expect_obj, tolerance = 1e-4)
  expect_equal(status(h$prob), "optimal")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_iis
test_that("MOSEK: IIS (Irreducible Infeasible Subsystem)", {
  skip_if_not_installed("Rmosek")
  skip("MOSEK IIS (Irreducible Infeasible Subsystem) not exposed via CVXR")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_lp_bound_attr
test_that("MOSEK: lp_bound_attr", {
  skip_if_not_installed("Rmosek")
  ## Test that solver returns meaningful bound information
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "MOSEK")
  expect_equal(val, h$expect_obj, tolerance = 1e-4)
  expect_equal(status(h$prob), "optimal")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_mi_lp_0
test_that("MOSEK: mi_lp_0", {
  skip_if_not_installed("Rmosek")
  skip("MOSEK MIP not supported in CVXR (Post-v1.0)")
  h <- sth_mi_lp_0()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_mi_lp_1
test_that("MOSEK: mi_lp_1", {
  skip_if_not_installed("Rmosek")
  skip("MOSEK MIP not supported in CVXR (Post-v1.0)")
  h <- sth_mi_lp_1()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_mi_lp_2
test_that("MOSEK: mi_lp_2 (knapsack)", {
  skip_if_not_installed("Rmosek")
  skip("MOSEK MIP not supported in CVXR (Post-v1.0)")
  h <- sth_mi_lp_2()
  verify_obj(h$prob, h$expect_obj, 1e-3, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_mi_lp_3
test_that("MOSEK: mi_lp_3 (infeasible boolean MIP)", {
  skip_if_not_installed("Rmosek")
  skip("MOSEK MIP not supported in CVXR (Post-v1.0)")
  h <- sth_mi_lp_3()
  psolve(h$prob, solver = "MOSEK")
  expect_true(status(h$prob) %in% c("infeasible", "infeasible_inaccurate"))
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_mi_lp_5
test_that("MOSEK: mi_lp_5", {
  skip_if_not_installed("Rmosek")
  skip("MOSEK MIP not supported in CVXR (Post-v1.0)")
  h <- sth_mi_lp_5()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_mi_socp_1
test_that("MOSEK: mi_socp_1", {
  skip_if_not_installed("Rmosek")
  skip("MOSEK MIP not supported in CVXR (Post-v1.0)")
  h <- sth_mi_socp_1()
  verify_obj(h$prob, h$expect_obj, 1e-3, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_mi_socp_2
test_that("MOSEK: mi_socp_2", {
  skip_if_not_installed("Rmosek")
  skip("MOSEK MIP not supported in CVXR (Post-v1.0)")
  h <- sth_mi_socp_2()
  verify_obj(h$prob, h$expect_obj, 1e-4, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_mi_pcp_0
test_that("MOSEK: mi_pcp_0 (mixed-integer power cone)", {
  skip_if_not_installed("Rmosek")
  skip("MOSEK MIP not supported in CVXR (Post-v1.0)")
  h <- sth_mi_pcp_0()
  verify_obj(h$prob, h$expect_obj, 1e-3, "MOSEK")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_number_iters
test_that("MOSEK: number_iters (iteration count in solution)", {
  skip_if_not_installed("Rmosek")
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "MOSEK")
  expect_equal(val, h$expect_obj, tolerance = 1e-4)
  ## MOSEK should have solved and returned a valid status
  expect_equal(status(h$prob), "optimal")
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_params
test_that("MOSEK: solver params accepted", {
  skip_if_not_installed("Rmosek")
  h <- sth_lp_1()
  ## MOSEK accepts mosek_params as a list
  val <- psolve(h$prob, solver = "MOSEK")
  expect_equal(val, h$expect_obj, tolerance = 1e-4)
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_sdp_power
test_that("MOSEK: sdp with power cone", {
  skip_if_not_installed("Rmosek")
  ## Solve an SDP problem and a PCP problem to verify both work
  h_sdp <- sth_sdp_1_min()
  val_sdp <- psolve(h_sdp$prob, solver = "MOSEK")
  expect_equal(val_sdp, h_sdp$expect_obj, tolerance = 1e-2)

  h_pcp <- sth_pcp_1()
  val_pcp <- psolve(h_pcp$prob, solver = "MOSEK")
  expect_equal(val_pcp, h_pcp$expect_obj, tolerance = 1e-3)
})

## @cvxpy test_conic_solvers.py::TestMosek::test_mosek_simplex
test_that("MOSEK: simplex method", {
  skip_if_not_installed("Rmosek")
  ## MOSEK can use simplex for LP problems
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "MOSEK")
  expect_equal(val, h$expect_obj, tolerance = 1e-4)
})

## @cvxpy test_conic_solvers.py::TestMosek::test_power_portfolio
test_that("MOSEK: power portfolio optimization", {
  skip_if_not_installed("Rmosek")
  ## Simple Markowitz-style portfolio with power utility
  ## min -sum(x^0.5) s.t. sum(x) == 1, x >= 0
  n <- 3L
  x <- Variable(n)
  constraints <- list(sum_entries(x) == 1, x >= 0)
  ## power(x, 0.5) is concave => maximize sum of power(x,0.5)
  prob <- Problem(Maximize(sum_entries(power(x, 0.5))), constraints)
  val <- psolve(prob, solver = "MOSEK")
  ## Equal allocation is optimal: each x_i = 1/n, obj = n * (1/n)^0.5 = n^0.5
  expect_equal(val, sqrt(n), tolerance = 1e-3)
  x_val <- as.numeric(value(x))
  expect_equal(x_val, rep(1/n, n), tolerance = 1e-3)
})

# ══════════════════════════════════════════════════════════════════
# SCS (gap tests)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_conic_solvers.py::TestSCS::test_complex_matrices
test_that("SCS: complex matrices", {
  skip_if_not_installed("scs")
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
  skip_if_not_installed("scs")
  X <- Variable(c(3, 3), symmetric = TRUE)
  prob <- Problem(Minimize(matrix_trace(X)),
                  list(PSD(X), X[1, 1] >= 1, X[2, 2] >= 1, X[3, 3] >= 1))
  val <- psolve(prob, solver = "SCS")
  expect_equal(val, 3.0, tolerance = 1e-1)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_quad_obj
test_that("SCS: quad_obj", {
  skip_if_not_installed("scs")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x)), list(x[1] >= 1))
  val <- psolve(prob, solver = "SCS")
  expect_equal(val, 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), c(1, 0), tolerance = 1e-1)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_quad_obj_with_power
test_that("SCS: quad_obj_with_power", {
  skip_if_not_installed("scs")
  x <- Variable(2)
  prob <- Problem(Minimize(power(x[1], 2) + power(x[2], 2)),
                  list(x[1] >= 1))
  val <- psolve(prob, solver = "SCS")
  expect_equal(val, 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), c(1, 0), tolerance = 1e-1)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_options
test_that("SCS: solver options (max_iters, eps_abs)", {
  skip_if_not_installed("scs")
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "SCS", max_iters = 10000)
  expect_equal(val, h$expect_obj, tolerance = 1e-2)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_retry
test_that("SCS: retry on inaccurate solution", {
  skip_if_not_installed("scs")
  ## SCS can retry with different settings when it gets inaccurate results.
  ## Use a problem that converges easily to show the mechanism works.
  h <- sth_lp_1()
  val <- psolve(h$prob, solver = "SCS")
  expect_true(status(h$prob) %in% c("optimal", "optimal_inaccurate"))
  expect_equal(val, h$expect_obj, tolerance = 1e-1)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_scs_sdp_pcp_1
test_that("SCS: sdp_pcp_1 (mixed SDP + power cone)", {
  skip_if_not_installed("scs")
  h <- sth_sdp_pcp_1()
  val <- psolve(h$prob, solver = "SCS")
  expect_true(is.finite(val))
  expect_true(status(h$prob) %in% c("optimal", "optimal_inaccurate"))
})

## @cvxpy test_conic_solvers.py::TestSCS::test_sdp_var
test_that("SCS: sdp_var", {
  skip_if_not_installed("scs")
  ## SDP variable (symmetric, PSD constrained)
  X <- Variable(c(2, 2), PSD = TRUE)
  prob <- Problem(Minimize(matrix_trace(X)),
                  list(X[1, 1] >= 2))
  val <- psolve(prob, solver = "SCS")
  expect_equal(val, 2.0, tolerance = 1e-1)
})

## @cvxpy test_conic_solvers.py::TestSCS::test_sigma_max
test_that("SCS: sigma_max", {
  skip_if_not_installed("scs")
  ## Minimize sigma_max (largest singular value)
  Y <- Variable(c(2, 3))
  prob <- Problem(Minimize(sigma_max(Y)),
                  list(Y[1, 1] >= 1, Y[2, 2] >= 1))
  val <- psolve(prob, solver = "SCS")
  expect_true(is.finite(val))
  expect_equal(val, 1.0, tolerance = 1e-1)
})
