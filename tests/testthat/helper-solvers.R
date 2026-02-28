## Solver availability guard for tests.
## Maps solver constant names to R package names via CVXR:::.SOLVER_PACKAGES.
## Usage: require_solver("GLPK_MI")  -- skips test if Rglpk not installed
require_solver <- function(solver_name) {
  pkg <- CVXR:::.SOLVER_PACKAGES[[solver_name]]
  if (is.null(pkg)) testthat::skip(paste0("Unavailable solver: ", solver_name))
  testthat::skip_if_not_installed(pkg)
}
