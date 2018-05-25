#'
#' The GLPK class
#'
#' This class is an interface for Gnu Linear Programming Toolkit solver
#'
#' @seealso the \href{https://www.gnu.org/software/glpk/}{Gnu GLKP site}.
#' @name GLPK-class
#' @aliases GLPK
#' @rdname GLPK-class
setClass("GLPK", contains = "Solver")



#' @rdname GLPK-class
#' @export
GLPK <- function() {
  new("GLPK")
}

#' @param object,solver A \linkS4class{GLPK} object.
#' @describeIn GLPK GLPK can handle linear programs.
setMethod("lp_capable", "GLPK", function(solver) { TRUE })

#' @describeIn GLPK GLPK can handle second-order cone programs.
setMethod("socp_capable", "GLPK", function(solver) { FALSE })

#' @describeIn GLPK GLPK can handle semidefinite programs.
setMethod("sdp_capable", "GLPK", function(solver) { FALSE })

#' @describeIn GLPK GLPK cannot handle exponential cone programs.
setMethod("exp_capable", "GLPK", function(solver) { FALSE })

#' @describeIn GLPK GLPK cannot handle mixed-integer programs.
setMethod("mip_capable", "GLPK", function(solver) { TRUE })

#'
#' GLPK Status Map
#'
#' Map of GLPK status to CVXR status.
#'
#' @param solver A \linkS4class{GLPK} object.
#' @param status An exit code returned by GLPK.
#' @return A string indicating the status, either "optimal", "infeasible", "unbounded", "optimal_inaccurate", "infeasible_inaccurate", "unbounded_inaccurate", or "solver_error".
#' @docType methods
#' @rdname GLPK-status_map
setMethod("status_map", "GLPK", function(solver, status) {

    ##GLP_UNDEF <- 1  /* solution is undefined */
    ##GLP_FEAS <- 2  /* solution is feasible */
    ##GLP_INFEAS <- 3  /* solution is infeasible */
    ##GLP_NOFEAS <- 4  /* no feasible solution exists */
    ##GLP_OPT <- 5  /* solution is optimal */
    ##GLP_UNBND <- 6  /* solution is unbounded */

    if(status == 1) {
        UNDEFINED
    } else if(status == 2) {
        OPTIMAL_INACCURATE
    } else if(status == 3 || status == 4) {
        INFEASIBLE
    } else if(status == 5) {
        OPTIMAL
    } else if(status == 6) {
        UNBOUNDED
    } else stop("GLPK status unrecognized: ", status)
})
#
#' @describeIn GLPK The name of the solver.
#' @export
setMethod("name", "GLPK", function(object) { GLPK_NAME })

#' @describeIn GLPK Imports the Rgpkk library.
setMethod("import_solver", "GLPK", function(solver) {
    requireNamespace("Rglpk", quietly = TRUE) &&
        requireNamespace("slam", quietly = TRUE)
})

setMethod("split_constr", "GLPK", function(solver, constr_map) {
    list(eq_constr = c(constr_map[[EQ_MAP]], constr_map[[LEQ_MAP]]), ineq_constr = list(), nonlin_constr = list())
})

#' @param objective A list representing the canonicalized objective.
#' @param constraints A list of canonicalized constraints.
#' @param cached_data A list mapping solver name to cached problem data.
#' @param warm_start A logical value indicating whether the previous solver result should be used to warm start.
#' @param verbose A logical value indicating whether to print solver output.
#' @param ... Additional arguments to the solver.
#' @describeIn GLPK Call the solver on the canonicalized problem.
setMethod("Solver.solve", "GLPK", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {
    solver_opts <- list(...)
    if (verbose) {
        solver_opts$verbose <- verbose
    }
    solver_opts$canonicalize_status <- FALSE

    data <- Solver.get_problem_data(solver, objective, constraints, cached_data)

    c <- data[[C_KEY]]
    dims <- data[[DIMS]]
    nvar <- length(c)
    A <- data[[A_KEY]]
    b <- data[[B_KEY]]

    bounds <- list(lower = list(ind = seq_along(c), val = rep(-Inf, nvar)))
    types <- rep("C", nvar)
    bools <- data[[BOOL_IDX]]
    ints <- data[[INT_IDX]]
    if (length(bools) > 0) {
        types[bools] <- "B"
    }
    if (length(ints) > 0) {
        types[ints] <- "I"
    }

    results_dict <- Rglpk::Rglpk_solve_LP(obj = c,
                                          mat = slam::as.simple_triplet_matrix(A),
                                          dir = c(rep("==", dims[[EQ_DIM]]),
                                                  rep("<=", dims[[LEQ_DIM]])),
                                          rhs = b,
                                          bounds = bounds,
                                          types = types,
                                          control = solver_opts,
                                          max = FALSE)
    format_results(solver, results_dict, data, cached_data)
})

#' @param results_dict A list containing the solver output.
#' @param data A list containing information about the problem.
#' @describeIn GLPK Convert raw solver output into standard list of results.
setMethod("format_results", "GLPK", function(solver, results_dict, data, cached_data) {
    new_results <- list()
    new_results[[STATUS]] <- status_map(solver, results_dict[[STATUS]])
    new_results[[SOLVE_TIME]] <- NA
    new_results[[SETUP_TIME]] <- NA

    if (new_results[[STATUS]] %in% SOLUTION_PRESENT) {
        ## Get primal variable values
        new_results[[PRIMAL]] <- results_dict[["solution"]]
        ## Get objective value
        new_results[[VALUE]] <- results_dict[["optimum"]] + data[[OFFSET]]
    }
    new_results
})

