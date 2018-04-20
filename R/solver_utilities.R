####################
#                  #
# Solver utilities #
#                  #
####################
# solver_intf <- list(ECOS(), ECOS_BB(), CVXOPT(), GLPK(), GLPK_MI(), CBC(), SCS(), GUROBI(), Elemental(), MOSEK(), LS())
solver_intf <- list(ECOS(), ECOS_BB(), SCS(), MOSEK(), LPSOLVE(), GLPK())
SOLVERS <- solver_intf
names(SOLVERS) <- sapply(solver_intf, function(solver) { name(solver) })

#'
#' Installed Solvers
#'
#' @return The names of all the installed solvers.
#' @docType methods
#' @rdname installed_solvers
#' @export
installed_solvers <- function() {
  installed <- sapply(SOLVERS, is_installed)
  names(SOLVERS)[which(installed)]

}

