####################
#                  #
# Solver utilities #
#                  #
####################
# solver_conic_intf <- list(ECOS(), ECOS_BB(), CVXOPT(), GLPK(), XPRESS(), GLPK_MI(), CBC_CONIC(), SCS(), SuperSCS(), GUROBI_CONIC(), MOSEK(), CPLEX_CONIC())
# solver_qp_intf <- list(OSQP(), GUROBI_QP(), CPLEX_QP())

solver_conic_intf <- list(ECOS(), ECOS_BB(), CBC_CONIC(), CPLEX_CONIC(), CVXOPT(),
                          GLPK_MI(), GLPK(), SCS(), GUROBI_CONIC(), MOSEK(), CLARABEL())
solver_qp_intf <- list(OSQP(), GUROBI_QP(), CPLEX_QP())

SOLVER_MAP_CONIC <- solver_conic_intf
names(SOLVER_MAP_CONIC) <- sapply(solver_conic_intf, name)

SOLVER_MAP_QP <- solver_qp_intf
names(SOLVER_MAP_QP) <- sapply(solver_qp_intf, name)

CONIC_SOLVERS <- c(MOSEK_NAME, ECOS_NAME, SUPER_SCS_NAME, SCS_NAME,
                   CPLEX_NAME, GUROBI_NAME, GLPK_NAME, XPRESS_NAME,
                   GLPK_MI_NAME, CBC_NAME, CVXOPT_NAME, ECOS_BB_NAME, CLARABEL_NAME)
QP_SOLVERS <- c(OSQP_NAME, GUROBI_NAME, CPLEX_NAME)

## Global variable for changing behavior
.CVXR_options <- new.env(parent = emptyenv())
.CVXR_options$blacklisted_solvers  <- character(0)

## CVXR global cache for storing solver codes etc.
.CVXR_cache <- new.env(parent = emptyenv())

#'
#' List installed solvers
#'
#' List available solvers, taking currently blacklisted solvers into
#' account.
#'
#' @return The names of all the installed solvers as a character vector.
#' @export
installed_solvers <- function() {
    ## Check conic solvers.
    installed_conic <- names(SOLVER_MAP_CONIC[sapply(SOLVER_MAP_CONIC, is_installed)])

    ## Check QP solvers.
    installed_qp <- names(SOLVER_MAP_QP[sapply(SOLVER_MAP_QP, is_installed)])

    setdiff(c(installed_conic, installed_qp), .CVXR_options$blacklisted_solvers)
}

#'
#' @param solvers a character vector of solver names, default \code{character(0)}
#' @return The current blacklist (character vector), invisibly.
#' @describeIn installed_solvers Add to solver blacklist, useful for temporarily disabling a solver
#' @export
add_to_solver_blacklist <- function(solvers) {
    stopifnot(is.character(solvers))
    result <- unique(c(.CVXR_options$blacklisted_solvers, solvers))
    .CVXR_options$blacklisted_solvers  <- result
    invisible(result)
}

#'
#' @describeIn installed_solvers Remove solvers from blacklist
#' @export
remove_from_solver_blacklist <- function(solvers) {
    stopifnot(is.character(solvers))
    result <- setdiff(.CVXR_options$blacklisted_solvers, solvers)
    .CVXR_options$blacklisted_solvers  <- result
    invisible(result)
}

#'
#' @describeIn installed_solvers Set solver blacklist to a value
#' @export
set_solver_blacklist <- function(solvers) {
    stopifnot(is.character(solvers))
    .CVXR_options$blacklisted_solvers  <- solvers
    invisible(solvers)
}

###
### Get solver codes from a saved file in `extdata` using solver
### canonical name and cache it if not already. For internal use
### @param name canonical name for solver, see utilities.R for names of solvers
### @return data frame of at least four columns `status` (character),
###   `code` (int), `cvxr_status`, and `description` (character).
get_solver_codes <- function(name) {
  result <- .CVXR_cache$status_codes[[name]]
  if (is.null(result)) {
    result <- .CVXR_cache$status_codes[[name]] <- 
      utils::read.csv(system.file("extdata", paste0(name, "_status_codes.csv"), package = "CVXR"))
  }
  result
}
    
