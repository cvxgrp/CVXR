## Cite verification pass, side-finding group 2 (notes/session_handoff_2026-10-07_cvxr_1.9.3_post_gurobi.md):
## upstream message wording at psd.py:56 and eval_params.py:15, and the pre-solve
## parameter check of problem.py:1019-1024, CVXR 1.9.2.9175. Expected messages were
## produced by cvxpy 1.9.3 in a pinned `uv run` process (scratch g2b_oracle.py). Every
## block was seen failing on 1.9.2.9174.

## cli wraps long messages at the console width, so compare with whitespace collapsed.
.expect_unset_param <- function(expr, nm) {
  cond <- expect_error(expr, class = "ParameterError")
  expect_identical(
    gsub("\\s+", " ", conditionMessage(cond)),
    sprintf("A Parameter (whose name is '%s') does not have a value associated with it; all Parameter objects must have values before solving a problem.", nm))
}

## @cvxpy test_constraints.py::TestConstraints::test_psd_constraint
test_that("PSD: a non-square argument is rejected in upstream words", {
  expect_error(PSD(Variable(c(2L, 3L))),
               "Non-square matrix in positive definite constraint.", fixed = TRUE)
})

## @cvxpy test_problem.py::TestProblem::test_solving_a_problem_with_unspecified_parameters
test_that("psolve: an unset Parameter is rejected before any solve path, naming it", {
  require_solver("CLARABEL")
  require_solver("SCS")
  q <- Parameter(name = "q"); x <- Variable()
  prob <- Problem(Minimize(q * x), list(x >= 1))
  .expect_unset_param(psolve(prob, solver = "CLARABEL"), "q")
  ## solver_path: ParameterError is not a SolverError, so it escapes the loop
  ## (problem.py:587 catches only SolverError).
  .expect_unset_param(psolve(prob, solver_path = c("CLARABEL", "SCS")), "q")
})

## @cvxpy test_problem.py::TestProblem::test_solving_a_problem_with_unspecified_parameters
test_that("psolve: a problem with no variables and an unset Parameter raises", {
  ## 1.9.2.9174 returned NULL with status "optimal" here.
  p <- Parameter(name = "lambda")
  .expect_unset_param(psolve(Problem(Minimize(p), list())), "lambda")
})

## @cvxpy test_dgp_dpp.py::TestDgpVariableBounds::test_solve_without_param_value_raises
test_that("psolve: the check runs on the gp and qcp paths too", {
  require_solver("CLARABEL")
  pp <- Parameter(pos = TRUE, name = "pp"); y <- Variable(pos = TRUE)
  .expect_unset_param(
    psolve(Problem(Minimize(pp * y), list(y >= 2)), solver = "CLARABEL", gp = TRUE), "pp")
  z <- Variable(); r <- Parameter(name = "r")
  .expect_unset_param(
    psolve(Problem(Minimize(ceil_expr(z)), list(z >= r)), solver = "CLARABEL", qcp = TRUE), "r")
})

## @cvxpy NONE
test_that("problem_data: the non-DPP path raises in eval_params.py:15 words", {
  require_solver("CLARABEL")
  q1 <- Parameter(name = "q1"); q2 <- Parameter(name = "q2"); t1 <- Variable()
  prob <- Problem(Minimize(q1 * q2 * t1), list(t1 >= 1))
  expect_error(suppressWarnings(problem_data(prob, solver = "CLARABEL")),
               "Problem contains unspecified parameters.", fixed = TRUE,
               class = "ParameterError")
})
