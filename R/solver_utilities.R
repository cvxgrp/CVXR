####################
#                  #
# Solver utilities #
#                  #
####################
solver_conic_intf <- list(ECOS(), ECOS_BB(), CVXOPT(), GLPK(), XPRESS(), GLPK_MI(), CBC_CONIC(), SCS(), SuperSCS(), GUROBI_CONIC(), MOSEK(), CPLEX_CONIC())
solver_qp_intf <- list(OSQP(), GUROBI_QP(), CPLEX_QP())

SOLVER_MAP_CONIC <- solver_conic_intf
names(SOLVER_MAP_CONIC) <- sapply(solver_conic_intf, name)

SOLVER_MAP_QP <- solver_qp_intf
names(SOLVER_MAP_QP) <- sapply(solver_qp_intf, name)

CONIC_SOLVERS <- c(MOSEK_NAME, ECOS_NAME, SUPER_SCS_NAME, SCS_NAME, 
                   CPLEX_NAME, GUROBI_NAME, GLPK_NAME, XPRESS_NAME, 
                   GLPK_MI_NAME, CBC_NAME, CVXOPT_NAME, ECOS_BB_NAME)
QP_SOLVERS <- c(OSQP_NAME, GUROBI_NAME, CPLEX_NAME)

#'
#' Installed Solvers
#'
#' @return The names of all the installed solvers.
#' @docType methods
#' @rdname installed_solvers
#' @export
installed_solvers <- function() {
  installed <- c()
  
  # Check conic solvers.
  installed_conic <- sapply(SOLVER_MAP_CONIC, is_installed)
  installed <- c(installed, names(SOLVER_MAP_CONIC)[installed_conic])
  
  # Check QP solvers.
  installed_qp <- sapply(SOLVER_MAP_QP, is_installed)
  installed <- c(installed, names(SOLVER_MAP_QP)[installed_qp])
  
  return(unique(installed))
}

INSTALLED_SOLVERS <- installed_solvers()
INSTALLED_CONIC_SOLVERS <- INSTALLED_SOLVERS[INSTALLED_SOLVERS %in% CONIC_SOLVERS]
