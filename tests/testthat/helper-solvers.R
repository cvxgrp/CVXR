## Solver availability guard for tests.
## Delegates to CVXR:::.solver_package_available() so the package owns the
## "is this solver usable?" policy (including the Rmosek >= 10 stub guard).
## Usage: require_solver("GLPK_MI")  -- skips test if Rglpk not installed
require_solver <- function(solver_name) {
  if (!CVXR:::.solver_package_available(solver_name)) {
    testthat::skip(paste0("Solver unavailable: ", solver_name))
  }
}
