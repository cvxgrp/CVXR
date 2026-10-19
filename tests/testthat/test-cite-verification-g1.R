## Cite verification pass, side-finding group 1, last item
## (notes/session_handoff_2026-10-07_cvxr_1.9.3_post_gurobi.md): with no usable solver,
## solving_chain.R raises upstream's "No installed solver could handle this problem."
## (solving_chain.py:101-103), CVXR 1.9.2.9181. The expected message was produced by
## cvxpy 1.9.3 in a pinned `uv run` process with every solver patched out
## (notes/audit/repro/39_no_installed_solver_message.py). Seen failing on 1.9.2.9180.

## @cvxpy NONE
test_that("auto-select with every solver excluded raises upstream's no-solver message", {
  old <- CVXR:::.cvxr_env$excluded_solvers
  on.exit(set_excluded_solvers(old), add = TRUE)
  set_excluded_solvers(installed_solvers())
  x <- Variable()
  expect_error(psolve(Problem(Minimize(x), list(x >= 1))),
               "No installed solver could handle this problem.",
               fixed = TRUE, class = "SolverError")
  expect_error(psolve(Problem(Maximize(log(x)), list(x <= 1))),
               "No installed solver could handle this problem.",
               fixed = TRUE, class = "SolverError")
})

## Mixed-integer problem, no solver named, no installed solver can take it: upstream's
## _fallback_solver (solving_chain.py:84-100) warns when the problem is not an LP, then
## raises "You need a mixed-integer solver for this model." CVXR 1.9.2.9192; the five
## cases were measured on cvxpy 1.9.3 first (scratch g9m/mi_probe.py: all five raise
## that error, C3 the MILP without the warning). On 1.9.2.9191 they gave the cone
## message, the MIQP message, or "could handle".
.mi_cases <- function() {
  y <- Variable(2, integer = TRUE)
  list(
    misocp = Problem(Maximize(sum(y)), list(p_norm(y, 2) <= 3)),
    milp   = Problem(Maximize(sum(y)), list(y <= 2.5)),
    miqp   = Problem(Minimize(sum_squares(y - 0.4)), list(y <= 2.5))
  )
}
.mi_message <- "You need a mixed-integer solver for this model."
.mi_warning <- "Your problem is mixed-integer but not an LP."

## @cvxpy NONE
test_that("mixed-integer, no capable solver installed (OSQP+HIGHS+CLARABEL+SCS): upstream's message", {
  rlang::local_options(rlib_warning_verbosity = "verbose")
  old <- CVXR:::.cvxr_env$excluded_solvers
  on.exit(set_excluded_solvers(old), add = TRUE)
  set_excluded_solvers(setdiff(installed_solvers(), c("OSQP", "HIGHS", "CLARABEL", "SCS")))
  cs <- .mi_cases()
  for (nm in c("misocp", "miqp")) {
    expect_warning(
      expect_error(psolve(cs[[nm]]), .mi_message, fixed = TRUE, class = "SolverError"),
      .mi_warning, fixed = TRUE)
  }
})

## @cvxpy NONE
test_that("mixed-integer, no capable solver installed (OSQP+CLARABEL+SCS): upstream's message", {
  rlang::local_options(rlib_warning_verbosity = "verbose")
  old <- CVXR:::.cvxr_env$excluded_solvers
  on.exit(set_excluded_solvers(old), add = TRUE)
  set_excluded_solvers(setdiff(installed_solvers(), c("OSQP", "CLARABEL", "SCS")))
  cs <- .mi_cases()
  expect_no_warning(
    expect_error(psolve(cs$milp), .mi_message, fixed = TRUE, class = "SolverError"))
  for (nm in c("misocp", "miqp")) {
    expect_warning(
      expect_error(psolve(cs[[nm]]), .mi_message, fixed = TRUE, class = "SolverError"),
      .mi_warning, fixed = TRUE)
  }
})
