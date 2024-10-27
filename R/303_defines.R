## CVXPY SOURCE: cvxpy/reductions/solvers/defines.py

solver_conic_intf <-
  list(
    #, DIFFCP()
    ECOS()
  , CVXOPT()
  , GLPK()
    #, COPT(),
  , GLPK_MI()
  , GUROBI_CONIC()
    # , CBC()
  , CLARABEL()
  , SCS()
    #, SDPA(),
    #, GUROBI(),
  , MOSEK()
  , CPLEX_CONIC()
    #, NAG(),
    #, XPRESS(),
    #, SCIP(),
    #, SCIPY(),
    #, GLOP()
    #, PDLP(),
  , ECOS_BB()
)

solver_qp_intf <-
  list(
    OSQP()
  , GUROBI_QP()
  , CPLEX_QP()
    #, XPRESS()
    #, COPT(),
    #, PIQP()
    #, PROXQP()
  )

SOLVER_MAP_CONIC <- solver_conic_intf
names(SOLVER_MAP_CONIC) <- sapply(solver_conic_intf, name)

SOLVER_MAP_QP <- solver_qp_intf
names(SOLVER_MAP_QP) <- sapply(solver_qp_intf, name)

# CONIC_SOLVERS and QP_SOLVERS are sorted in order of decreasing solver
# preference. QP_SOLVERS are those for which we have written interfaces
# and are supported by QpSolver.
CONIC_SOLVERS <-
  c(
    MOSEK_NAME
  , ECOS_NAME
  , CLARABEL_NAME
  , SCS_NAME
    #, SDPA_NAME
  , CPLEX_NAME
  , GUROBI_NAME
    #, COPT_NAME
  , GLPK_NAME
    #, NAG_NAME
  , GLPK_MI_NAME
  , CBC_NAME
  , CVXOPT_NAME
    #, XPRESS_NAME
    #, DIFFCP_NAME
    #, SCIP_NAME
    #, GLOP_NAME
    #, PDLP_NAME
  , ECOS_BB_NAME
  )

QP_SOLVERS <-
  c(
    OSQP_NAME
  , GUROBI_NAME
  , CPLEX_NAME
    #, PIQP_NAME
    #, XPRESS_NAME
    #, COPT_NAME
    #, PROXQP_NAME
  )

MI_SOLVERS <-
  c(
    GLPK_MI_NAME
  , MOSEK_NAME
  , GUROBI_NAME
  , CPLEX_NAME
    #, XPRESS_NAME
  , CBC_NAME
    #, SCIP_NAME
    #, COPT_NAME
  , ECOS_BB_NAME
)

MI_SOCP_SOLVERS <-
  c(
    MOSEK_NAME
  , GUROBI_NAME
  , CPLEX_NAME
    #, XPRESS_NAME
    #, SCIP_NAME
  , ECOS_BB_NAME
  )

## Global variable for changing behavior
#.CVXR_options <- new.env(parent = emptyenv())
.CVXR_options$blacklisted_solvers  <- character(0)

#' @param solvers a character vector of solver names
#' @describeIn installed_solvers Add to solver blacklist, useful for temporarily disabling a solver
#' @export
add_to_solver_blacklist <- function(solvers) {
    stopifnot(is.character(solvers))
    result <- unique(c(.CVXR_options$blacklisted_solvers, solvers))
    .CVXR_options$blacklisted_solvers  <- result
    invisible(result)
}

#' @describeIn installed_solvers Remove solvers from blacklist
#' @export
remove_from_solver_blacklist <- function(solvers) {
    stopifnot(is.character(solvers))
    result <- setdiff(.CVXR_options$blacklisted_solvers, solvers)
    .CVXR_options$blacklisted_solvers  <- result
    invisible(result)
}

#' @describeIn installed_solvers Set solver blacklist to a value
#' @export
set_solver_blacklist <- function(solvers) {
    stopifnot(is.character(solvers))
    .CVXR_options$blacklisted_solvers  <- solvers
    invisible(solvers)
}

#'
#' List installed solvers, taking blacklisted solvers into account
#' @return The names of all the installed solvers as a character vector
#' @export
installed_solvers <- function() {
    ## Check conic solvers.
    installed_conic <- names(SOLVER_MAP_CONIC[sapply(SOLVER_MAP_CONIC, is_installed)])

    ## Check QP solvers.
    installed_qp <- names(SOLVER_MAP_QP[sapply(SOLVER_MAP_QP, is_installed)])

    # Remove duplicate names (for solvers that handle both conic and QP)
    installed <- unique(c(installed_conic, installed_qp))

    # Remote blacklisted solvers
    setdiff(installed, .CVXR_options$blacklisted_solvers)
}

INSTALLED_SOLVERS <- installed_solvers()
INSTALLED_CONIC_SOLVERS <- intersect(INSTALLED_SOLVERS, CONIC_SOLVERS)
INSTALLED_MI_SOLVERS <- intersect(INSTALLED_SOLVERS, MI_SOLVERS)

