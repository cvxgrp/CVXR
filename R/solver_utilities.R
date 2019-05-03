####################
#                  #
# Solver utilities #
#                  #
####################
# solver_intf <- list(ECOS(), ECOS_BB(), CVXOPT(), GLPK(), GLPK_MI(), CBC(), SCS(), GUROBI(), Elemental(), MOSEK(), LS())
solver_intf <- list(ECOS(), ECOS_BB(), SCS(), MOSEK(), LPSOLVE(), GLPK(), GUROBI_CONIC())
solver_conic_intf <- list(ECOS(), ECOS_BB(), CVXOPT(), GLPK(), XPRESS(), GLPK_MI(), CBC_CONIC(), SCS(), SuperSCS(), GUROBI_CONIC(), Elemental(), MOSEK(), CPLEX_CONIC())
solver_qp_intf <- list(OSQP(), GUROBI_QP(), CPLEX_QP())

SOLVERS <- solver_intf
names(SOLVERS) <- sapply(solver_intf, name)

SOLVER_MAP_CONIC <- solver_conic_intf
names(SOLVER_MAP_CONIC) <- sapply(solver_conic_intf, name)

SOLVER_MAP_QP <- solver_qp_intf
names(SOLVER_MAP_QP) <- sapply(solver_qp_intf, name)

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

