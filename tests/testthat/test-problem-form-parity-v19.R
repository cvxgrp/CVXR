## CVXPY 1.9 ProblemForm / solving-chain parity.
##
## CVXR does not expose CVXPY's ProblemForm class directly. These tests cover
## the equivalent solver-chain predicates, cone detection/stuffing, and GP
## dispatch behavior that CVXPY's test_problem_form.py exercises.

problem_form_cone_names <- function(prob) {
  vapply(CVXR:::.required_cone_types(prob), function(cls) cls@name, character(1L))
}

problem_form_chain_solver <- function(prob, solver = NULL, gp = FALSE) {
  chain <- CVXR:::construct_solving_chain(prob, solver = solver, gp = gp)
  CVXR:::solver_name(chain@solver)
}

problem_form_is_commercial_solver <- function(solver) {
  solver %in% c("MOSEK", "GUROBI", "CPLEX", "XPRESS")
}

problem_form_expect_solver <- function(actual, expected) {
  expect_true(actual == expected || problem_form_is_commercial_solver(actual),
              info = sprintf("expected %s or commercial solver, got %s",
                             expected, actual))
}

problem_form_is_quad_stuffing <- function(reduction) {
  S7::S7_inherits(reduction, ConeMatrixStuffing) && isTRUE(reduction@quad_obj)
}

## @cvxpy test_problem_form.py::TestProblemForm::test_lp
test_that("ProblemForm parity: LP structure has linear objective and base cones", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1, x[1] == 2))

  expect_false(CVXR:::has_quadratic_term(prob@objective@args[[1L]]))
  expect_false(is_mixed_integer(prob))
  expect_true(length(prob@constraints) > 0L)

  data <- problem_data(prob, solver = "CLARABEL")$data
  dims <- data[[CVXR:::SD_DIMS]]
  expect_equal(dims@zero, 1L)
  expect_equal(dims@nonneg, 2L)
  expect_length(problem_form_cone_names(prob), 0L)
})

## @cvxpy test_problem_form.py::TestProblemForm::test_qp
test_that("ProblemForm parity: QP structure detects quadratic objective", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x)), list(x >= 1))

  expect_true(CVXR:::has_quadratic_term(prob@objective@args[[1L]]))
  expect_true(length(prob@constraints) > 0L)
})

## @cvxpy test_problem_form.py::TestProblemForm::test_qp_filtering
test_that("ProblemForm parity: QP path filters quadratic SOCs but keeps other cones", {
  x <- Variable()
  qp <- Problem(Minimize(sum_squares(x)))
  qp_cone <- reduction_apply(CVXR:::Dcp2Cone(quad_obj = TRUE), qp)[[1L]]
  expect_false(any(vapply(qp_cone@constraints, S7::S7_inherits, logical(1L), SOC)))

  y <- Variable()
  exp_qp <- Problem(Minimize(quad_over_lin(exp(y), 1)))
  exp_cone <- reduction_apply(CVXR:::Dcp2Cone(quad_obj = TRUE), exp_qp)[[1L]]
  expect_true(any(vapply(exp_cone@constraints, S7::S7_inherits, logical(1L), ExpCone)))
  expect_false(any(vapply(exp_cone@constraints, S7::S7_inherits, logical(1L), SOC)))

  z <- Variable()
  nested_qp <- Problem(Minimize(sum_squares(sum_squares(z))))
  nested_cone <- reduction_apply(CVXR:::Dcp2Cone(quad_obj = TRUE), nested_qp)[[1L]]
  expect_true(any(vapply(nested_cone@constraints, S7::S7_inherits, logical(1L), SOC)))
})

## @cvxpy test_problem_form.py::TestProblemForm::test_quad_over_lin_constant_denom_qp_path
test_that("ProblemForm parity: quad_over_lin with constant denominator stays on QP path", {
  x <- Variable()
  prob <- Problem(Minimize(quad_over_lin(exp(x), 1)))
  conic_prob <- reduction_apply(CVXR:::Dcp2Cone(quad_obj = TRUE), prob)[[1L]]

  expect_true(any(vapply(conic_prob@constraints, S7::S7_inherits, logical(1L), ExpCone)))
  expect_false(any(vapply(conic_prob@constraints, S7::S7_inherits, logical(1L), SOC)))
})

## @cvxpy test_problem_form.py::TestProblemForm::test_socp
test_that("ProblemForm parity: SOCP structure reports SOC cones", {
  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 2)), list(x >= 1))

  expect_true("SOC" %in% problem_form_cone_names(prob))
  expect_true(length(prob@constraints) > 0L)
})

## @cvxpy test_problem_form.py::TestProblemForm::test_sdp
test_that("ProblemForm parity: SDP structure reports PSD cones", {
  X <- Variable(c(2, 2), symmetric = TRUE)
  prob <- Problem(Minimize(lambda_max(X)), list(X %>>% 0))

  expect_true("PSD" %in% problem_form_cone_names(prob))
  expect_true(length(prob@constraints) > 0L)
})

## @cvxpy test_problem_form.py::TestProblemForm::test_exp_cone
test_that("ProblemForm parity: exponential atoms report ExpCone", {
  x <- Variable(2, pos = TRUE)
  prob <- Problem(Maximize(sum_entries(log(x))), list(x <= 2))

  expect_true("ExpCone" %in% problem_form_cone_names(prob))
})

## @cvxpy test_problem_form.py::TestProblemForm::test_pow_cone
test_that("ProblemForm parity: explicit power-cone constraints report PowCone3D", {
  x <- Variable()
  y <- Variable()
  z <- Variable()
  prob <- Problem(Minimize(x), list(PowCone3D(x, y, z, 0.5)))

  expect_true("PowCone3D" %in% problem_form_cone_names(prob))
})

## @cvxpy test_problem_form.py::TestProblemForm::test_mixed_integer
test_that("ProblemForm parity: integer and boolean variables are mixed-integer", {
  x <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 0, x <= 5))
  expect_true(is_mixed_integer(prob))

  y <- Variable(2, boolean = TRUE)
  prob2 <- Problem(Minimize(sum_entries(y)), list(y >= 0))
  expect_true(is_mixed_integer(prob2))
})

## @cvxpy test_problem_form.py::TestProblemForm::test_unconstrained_linear
test_that("ProblemForm parity: unconstrained linear objective has no cones", {
  x <- Variable()
  prob <- Problem(Minimize(x))

  expect_false(CVXR:::has_quadratic_term(prob@objective@args[[1L]]))
  expect_length(problem_form_cone_names(prob), 0L)
  expect_false(is_mixed_integer(prob))
  expect_length(prob@constraints, 0L)
})

## @cvxpy test_problem_form.py::TestProblemForm::test_sdp_nsd_variable
test_that("ProblemForm parity: NSD variables create PSD cone dimensions", {
  X <- Variable(c(2, 2), NSD = TRUE)
  prob <- Problem(Maximize(matrix_trace(X)))
  data <- problem_data(prob, solver = "SCS")$data

  expect_equal(data[[CVXR:::SD_DIMS]]@psd, 2L)
})

## @cvxpy test_problem_form.py::TestPickDefaultSolver::test_lp_gets_clarabel test_problem_form.py::TestPickDefaultSolver::test_qp_gets_osqp test_problem_form.py::TestPickDefaultSolver::test_sdp_gets_scs test_problem_form.py::TestPickDefaultSolver::test_socp_gets_clarabel test_problem_form.py::TestPickDefaultSolver::test_mi_gets_highs_or_premium
test_that("ProblemForm parity: default solver selection matches CVXPY policy modulo commercial solvers", {
  x <- Variable(2)
  lp <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  problem_form_expect_solver(problem_form_chain_solver(lp), "CLARABEL")

  q <- Variable(2)
  qp <- Problem(Minimize(sum_squares(q)), list(q >= 1))
  problem_form_expect_solver(problem_form_chain_solver(qp), "OSQP")

  X <- Variable(c(2, 2), symmetric = TRUE)
  sdp <- Problem(Minimize(lambda_max(X)), list(X %>>% 0))
  problem_form_expect_solver(problem_form_chain_solver(sdp), "SCS")

  s <- Variable(2)
  socp <- Problem(Minimize(p_norm(s, 2)), list(s >= 1))
  problem_form_expect_solver(problem_form_chain_solver(socp), "CLARABEL")

  m <- Variable(2, integer = TRUE)
  mip <- Problem(Minimize(sum_entries(m)), list(m >= 0, m <= 5))
  mi_solver <- problem_form_chain_solver(mip)
  expect_true(mi_solver == "HIGHS" || problem_form_is_commercial_solver(mi_solver),
              info = sprintf("expected HIGHS or commercial solver, got %s", mi_solver))
})

## @cvxpy test_problem_form.py::TestResolveAndBuildChain::test_solver_none_picks_default
test_that("ProblemForm parity: solver NULL builds a chain and solves LP", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  val <- psolve(prob, solver = NULL)
  expect_equal(val, 2, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-5)
})

## @cvxpy test_problem_form.py::TestResolveAndBuildChain::test_solver_string_clarabel
test_that("ProblemForm parity: named CLARABEL chain solves SOCP", {
  skip_if_not("CLARABEL" %in% installed_solvers())

  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 2)), list(x >= 1))
  expect_equal(problem_form_chain_solver(prob, solver = "CLARABEL"), "CLARABEL")

  val <- psolve(prob, solver = "CLARABEL")
  expect_equal(val, sqrt(2), tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-5)
})

## @cvxpy test_problem_form.py::TestResolveAndBuildChain::test_solver_not_installed
test_that("ProblemForm parity: unknown solver is rejected", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  expect_error(CVXR:::construct_solving_chain(prob, solver = "NONEXISTENT"),
               class = "rlang_error")
})

## @cvxpy test_problem_form.py::TestResolveAndBuildChain::test_solver_cannot_handle
test_that("ProblemForm parity: incompatible solver is rejected", {
  skip_if_not("OSQP" %in% installed_solvers())

  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 2)), list(x >= 1))

  expect_error(CVXR:::construct_solving_chain(prob, solver = "OSQP"),
               "does not support|cannot handle|required cone")
})

## @cvxpy test_problem_form.py::TestGP::test_gp_cones
test_that("ProblemForm parity: GP chain inserts Dgp2Dcp and preserves conic routing", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x + y), list(x * y >= 1))
  chain <- CVXR:::construct_solving_chain(prob, gp = TRUE)

  expect_true(any(vapply(chain@reductions, S7::S7_inherits, logical(1L), Dgp2Dcp)))
  expect_false(any(vapply(chain@reductions, problem_form_is_quad_stuffing, logical(1L))))

  x2 <- Variable(2, pos = TRUE)
  geo_prob <- Problem(Maximize(geo_mean(x2)), list(x2 <= 2))
  geo_chain <- CVXR:::construct_solving_chain(geo_prob, gp = TRUE)
  expect_true(any(vapply(geo_chain@reductions, S7::S7_inherits, logical(1L), Dgp2Dcp)))
})

## @cvxpy test_problem_form.py::TestGP::test_gp_solve_default
test_that("ProblemForm parity: GP with default solver solves correctly", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x + y), list(x * y >= 1))

  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 2, tolerance = 1e-3)
})

## @cvxpy test_problem_form.py::TestGP::test_gp_solve_named
test_that("ProblemForm parity: GP with named SCS solver solves correctly", {
  skip_if_not("SCS" %in% installed_solvers())

  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x + y), list(x * y >= 1))

  expect_equal(problem_form_chain_solver(prob, solver = "SCS", gp = TRUE), "SCS")
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(val, 2, tolerance = 1e-2)
})

# ===================================================================
# Solver resolution (ADR D_19.8). Each test is one class of the selection
# matrix, notes/solver_selection_matrix_2026-10-14.md, measured against CVXPY
# 1.9.3 (docsreg_9204/selmatrix.{py,R}).
# ===================================================================

## Run `code` with exactly `excluded` excluded and the commercial preference
## set to `commercial`, restoring both afterwards.
problem_form_with_solvers <- function(excluded, commercial, code) {
  old_ex <- CVXR:::.cvxr_env$excluded_solvers
  old_opt <- options(CVXR.default_to_commercial_solvers = commercial)
  on.exit({ set_excluded_solvers(if (is.null(old_ex)) character(0) else old_ex)
            options(old_opt) }, add = TRUE)
  set_excluded_solvers(excluded)
  for (id in c("cvxr_mixed_integer_not_lp",
               paste0("cvxr_default_solvers_unavailable_", installed_solvers())))
    rlang::reset_warning_verbosity(id)
  force(code)
}

problem_form_iface <- function(prob, solver = NULL) {
  sv <- problem_data(prob, solver = solver)$chain@solver
  if (S7::S7_inherits(sv, CVXR:::QpSolver)) "qp" else "conic"
}

## Matrix D3 (solving_chain.py:69-83): when no default solver is available,
## the fallback warns and names the solver it took. Upstream mocks
## pick_default_solver to return None; here the default solvers are excluded.
## @cvxpy test_problem_form.py::TestResolveAndBuildChain::test_fallback_when_default_missing
test_that("the fallback solver warns that the default solvers are unavailable", {
  require_solver("ECOS")
  x <- Variable(2)
  prob <- Problem(Minimize(sum(x)), list(x >= 1))
  problem_form_with_solvers(
    c("CLARABEL", "OSQP", "SCS", "HIGHS", "SCIP"), FALSE, {
      expect_warning(val <- psolve(prob),
                     "default solvers are not available.*ECOS")
      expect_equal(solver_stats(prob)@solver_name, "ECOS")
      expect_equal(val, 2, tolerance = 1e-6)
      expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-6)
    })
})

## Matrix D4 (solving_chain.py:71-76): the fallback walks INSTALLED_SOLVERS,
## conic interface before QP for each name, so a QP lands on ECOS, not PIQP.
## @cvxpy NONE
test_that("the fallback solver tries each name's conic interface first", {
  require_solver("ECOS")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - 1)), list(x >= 0))
  problem_form_with_solvers(
    c("CLARABEL", "OSQP", "SCS", "HIGHS", "SCIP", "MOSEK", "GUROBI", "CPLEX",
      "XPRESS"), TRUE, {
      expect_warning(val <- psolve(prob), "default solvers are not available")
      expect_equal(solver_stats(prob)@solver_name, "ECOS")
      expect_equal(val, 0, tolerance = 1e-6)
    })
})

## Matrix D1 (problem_form.py:356-360): the commercial-first preference covers
## QPs and MIQPs too, not only conic problems.
## @cvxpy NONE
test_that("the commercial-first default applies to a QP and an MIQP", {
  require_solver("MOSEK")
  x <- Variable(2)
  qp <- Problem(Minimize(sum_squares(x - 1)), list(x >= 0))
  y <- Variable(2, integer = TRUE)
  miqp <- Problem(Minimize(sum_squares(y - 0.4)), list(sum(y) >= 1))
  problem_form_with_solvers(character(0), TRUE, {
    expect_equal(psolve(qp), 0, tolerance = 1e-6)
    expect_equal(solver_stats(qp)@solver_name, "MOSEK")
    expect_equal(psolve(miqp), 0.52, tolerance = 1e-6)
    expect_equal(solver_stats(miqp)@solver_name, "MOSEK")
  })
})

## Matrix D2 and D7 (problem_form.py:362-371; no early MIQP rejection
## upstream): an MIQP goes to SCIP's conic interface, by default when no
## commercial solver is wanted, and when SCIP is named.
## @cvxpy NONE
test_that("an MIQP is solved by SCIP, by default and when named", {
  require_solver("SCIP")
  y <- Variable(2, integer = TRUE)
  miqp <- Problem(Minimize(sum_squares(y - 0.4)), list(sum(y) >= 1))
  problem_form_with_solvers(character(0), FALSE, {
    expect_equal(psolve(miqp), 0.52, tolerance = 1e-6)
    expect_equal(solver_stats(miqp)@solver_name, "SCIP")
  })
  y2 <- Variable(2, integer = TRUE)
  miqp2 <- Problem(Minimize(sum_squares(y2 - 0.4)), list(sum(y2) >= 1))
  expect_equal(psolve(miqp2, solver = "SCIP"), 0.52, tolerance = 1e-6)
})

## Matrix D6 (solving_chain.py:317-329): a named solver with both interfaces
## tries QP first only for a quadratic objective, so an LP goes conic. Duals
## measured equal to CVXPY 1.9.3 on all four (docsreg_9204/d6_duals.{R,py}).
## @cvxpy NONE
test_that("a named solver with both interfaces takes an LP through its conic interface", {
  for (sv in c("HIGHS", "GUROBI", "CPLEX", "XPRESS")) {
    if (!CVXR:::.solver_package_available(sv)) next
    x <- Variable(3)
    c2 <- sum(x) == 4
    c3 <- 2 * x[1] + x[2] <= 5
    prob <- Problem(Minimize(-3 * x[1] - x[2] + x[3]),
                    list(x >= c(1, 0, -1), c2, c3, x[3] <= 3))
    expect_equal(problem_form_iface(prob, sv), "conic", info = sv)
    expect_equal(psolve(prob, solver = sv), -6, tolerance = 1e-6, info = sv)
    expect_equal(as.numeric(dual_value(c2)), -1, tolerance = 1e-6, info = sv)
    expect_equal(as.numeric(dual_value(c3)), 2, tolerance = 1e-6, info = sv)
    q <- Variable(2)
    qp <- Problem(Minimize(sum_squares(q - 1)), list(q >= 0))
    expect_equal(problem_form_iface(qp, sv), "qp", info = sv)
  }
})

## RECORDED DIVERGENCE (user, 2026-10-14): solver selection honors
## use_quad_obj = FALSE, which CVXPY 1.9.3's can_solve ignores (it then picks
## a QP solver and fails on the objective's SOC; measured,
## notes/upstream_report_cvxpy_use_quad_obj.md). CVXR 1.9.2.9211 solved these.
## @cvxpy NONE
test_that("use_quad_obj = FALSE keeps a QP solver out of the selection", {
  x <- Variable(2)
  qp <- Problem(Minimize(sum_squares(x - 1)), list(x >= 0))
  problem_form_with_solvers(c("MOSEK", "GUROBI", "CPLEX", "XPRESS"), TRUE, {
    expect_equal(psolve(qp, use_quad_obj = FALSE), 0, tolerance = 1e-6)
    expect_equal(solver_stats(qp)@solver_name, "CLARABEL")
  })
  expect_error(psolve(qp, solver = "OSQP", use_quad_obj = FALSE),
               "does not support the required cone types: SOC",
               class = "SolverError")
  require_solver("GUROBI")
  x2 <- Variable(2)
  qp2 <- Problem(Minimize(sum_squares(x2 - 1)), list(x2 >= 0))
  expect_equal(psolve(qp2, solver = "GUROBI", use_quad_obj = FALSE), 0,
               tolerance = 1e-6)
  expect_false(S7::S7_inherits(
    problem_data(qp2, solver = "GUROBI", use_quad_obj = FALSE)$chain@solver,
    CVXR:::QpSolver))
})

## Matrix D8 (solver.py:168-169): REQUIRES_CONSTR is checked against the
## problem before canonicalization, so SCS refuses a problem with no
## constraints at selection time.
## @cvxpy NONE
test_that("SCS refuses a problem without constraints at solver selection", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum(x)))
  expect_error(psolve(prob, solver = "SCS"),
               "cannot solve this problem", class = "SolverError")
})

## MOSEK's MI_SUPPORTED_CONSTRAINTS (mosek_conif.py:114-115) has no PSD cone:
## MOSEK cannot solve a mixed-integer SDP, and says so at selection time.
## @cvxpy NONE
test_that("MOSEK is not selected for a mixed-integer SDP", {
  require_solver("MOSEK")
  X <- Variable(c(2, 2), PSD = TRUE)
  z <- Variable(integer = TRUE)
  prob <- Problem(Minimize(matrix_trace(X) + z), list(X[1, 1] >= 1, z >= 0.5))
  expect_error(psolve(prob, solver = "MOSEK"),
               "does not support the required cone types: PSD",
               class = "SolverError")
})

## solving_chain.py:71-73: the fallback never takes a commercial solver, so
## under the opt-out no commercial solver is chosen unless named. Measured vs
## cvxpy 1.9.3 (docsreg_9204/commercial_probe.{R,py}).
## @cvxpy NONE
test_that("with the commercial opt-out, auto-selection never picks a commercial solver", {
  require_solver("SCIP")
  y <- Variable(2, integer = TRUE)
  miqp <- Problem(Minimize(sum_squares(y - 0.4)), list(sum(y) >= 1))
  x <- Variable(integer = TRUE)
  mi_exp <- Problem(Minimize(exp(x)), list(x >= 1.5))
  problem_form_with_solvers(character(0), FALSE, {
    expect_equal(psolve(miqp), 0.52, tolerance = 1e-6)
    expect_false(solver_stats(miqp)@solver_name %in% CVXR:::.COMMERCIAL_CONIC_SOLVERS)
    ## SCIP has no exponential cone and ECOS_BB is never auto-selected.
    expect_error(suppressWarnings(psolve(mi_exp)),
                 "You need a mixed-integer solver", class = "SolverError")
  })
  ## A named commercial solver is unaffected by the opt-out.
  require_solver("MOSEK")
  x2 <- Variable(integer = TRUE)
  p2 <- Problem(Minimize(exp(x2)), list(x2 >= 1.5))
  problem_form_with_solvers(character(0), FALSE,
    expect_equal(psolve(p2, solver = "MOSEK"), exp(2), tolerance = 1e-6))
})

## RECORDED DIVERGENCE (user, 2026-10-14): ECOS_BB is never auto-selected.
## CVXPY's fallback can pick it, and on min exp(x), x >= 1.5, x integer it
## returns 1 with x = 0 and status "optimal" -- in CVXPY 1.9.3 (ecos 2.0.14)
## and in CVXR (ECOSolveR 0.6.1) alike (docsreg_9204/ecosbb_probe.{R,py}).
## Naming it still works.
## @cvxpy NONE
test_that("ECOS_BB is never auto-selected, only when named", {
  require_solver("ECOS_BB")
  x <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(sum(x)), list(x >= 1.5))
  problem_form_with_solvers(setdiff(installed_solvers(), c("ECOS", "ECOS_BB")),
                            TRUE, {
    expect_error(suppressWarnings(psolve(prob)),
                 "You need a mixed-integer solver", class = "SolverError")
    expect_equal(psolve(prob, solver = "ECOS_BB"), 4, tolerance = 1e-6)
    expect_equal(solver_stats(prob)@solver_name, "ECOS_BB")
  })
})
