## Solver availability guard for tests.
## Delegates to CVXR:::.solver_package_available() so the package owns the
## "is this solver usable?" policy (including the Rmosek >= 10 stub guard).
## Usage: require_solver("GLPK_MI")  -- skips test if Rglpk not installed
require_solver <- function(solver_name) {
  if (!CVXR:::.solver_package_available(solver_name)) {
    testthat::skip(paste0("Solver unavailable: ", solver_name))
  }
}

## Solvers CVXR does NOT Import. Imports: clarabel, highs, osqp, scs.
.NON_IMPORTS_SOLVERS <- c("MOSEK", "GUROBI", "CPLEX", "SCIP", "XPRESS",
                          "GLPK", "GLPK_MI", "ECOS", "ECOS_BB", "CVXOPT",
                          "PIQP")

## Run `code` with AUTOMATIC solver selection restricted to the Imports
## solvers, restoring the previous exclusion list afterwards.
##
## WHY, and why it is SCOPED rather than suite-wide. A maintainer's machine may
## have all 15 solvers licensed, so `psolve(prob)` with no `solver=` picks the
## highest-priority INSTALLED solver -- MOSEK here, which CI and CRAN never
## have. A test relying on auto-selection then measures a different solver for
## the maintainer than for anyone else. That is not theoretical: three DQCP
## bisection tests failed locally and passed on CI because auto-selection
## reached MOSEK and MOSEK could not solve the bisection subproblems. The
## failure looks exactly like a code regression, which is the expensive part.
##
## Excluding suite-wide is NOT the fix, and was tried: this package has whole
## files of tests that exist to exercise XPRESS, MIQP, MI-SOCP and other
## capabilities the Imports solvers do not have, and a blanket exclusion makes
## those error out with "no installed solver supports ...". Determinism is
## wanted for the tests that do not care which solver runs; the solver-specific
## tests must keep seeing their solver.
##
## So wrap only the tests whose POINT is the result, not the solver:
##     val <- with_imports_solvers(psolve(prob, qcp = TRUE))
##
## Exclusion governs AUTOMATIC selection only -- an explicit `solver = "ECOS"`
## inside the block still runs, and `installed_solvers()` still reports the
## full set, so `require_solver()` guards are unaffected.
with_imports_solvers <- function(code) {
  old <- CVXR:::.cvxr_env$excluded_solvers
  CVXR::set_excluded_solvers(.NON_IMPORTS_SOLVERS)
  ## on.exit, not a trailing call: the list must be restored even when `code`
  ## signals, or one failing test silently re-scopes every test after it.
  on.exit(CVXR::set_excluded_solvers(old), add = TRUE)
  force(code)
}
