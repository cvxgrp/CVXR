#'
#' The LPSOLVE class
#'
#' This class is an interface for Gnu Linear Programming Toolkit solver
#'
#' @seealso the \href{http://www.gnu.org/software/LPSOLVE/}{Gnu GLKP site}.
#' @name LPSOLVE-class
#' @aliases LPSOLVE
#' @rdname LPSOLVE-class
setClass("LPSOLVE", contains = "Solver")

#' @rdname LPSOLVE-class
#' @export
LPSOLVE <- function() {
  new("LPSOLVE")
}

#' @param object,solver A \linkS4class{LPSOLVE} object.
#' @describeIn LPSOLVE LPSOLVE can handle linear programs.
setMethod("lp_capable", "LPSOLVE", function(solver) { TRUE })

#' @describeIn LPSOLVE LPSOLVE can handle second-order cone programs.
setMethod("socp_capable", "LPSOLVE", function(solver) { FALSE })

#' @describeIn LPSOLVE LPSOLVE can handle semidefinite programs.
setMethod("sdp_capable", "LPSOLVE", function(solver) { FALSE })

#' @describeIn LPSOLVE LPSOLVE cannot handle exponential cone programs.
setMethod("exp_capable", "LPSOLVE", function(solver) { FALSE })

#' @describeIn LPSOLVE LPSOLVE cannot handle mixed-integer programs.
setMethod("mip_capable", "LPSOLVE", function(solver) { TRUE })

#'
#' LPSOLVE Status Map
#'
#' Map of LPSOLVE status to CVXR status.
#'
#' @param solver A \linkS4class{LPSOLVE} object.
#' @param status An exit code returned by LPSOLVE.
#' @return A string indicating the status, either "optimal", "infeasible", "unbounded", "optimal_inaccurate", "infeasible_inaccurate", "unbounded_inaccurate", or "solver_error".
#' @docType methods
#' @rdname LPSOLVE-status_map
setMethod("status_map", "LPSOLVE", function(solver, status) {
    if(status == 0) {
        OPTIMAL
    } else if(status == 1) {
        OPTIMAL_INACCURATE
    } else if(status == 2) {
        INFEASIBLE
    } else if(status == 3) {
        UNBOUNDED
    } else if (status %in% seq(4, 13)){
        SOLVER_ERROR
    } else {
        stop(sprintf("LPSOLVE status %d unrecognized!"))
    }
})
#
#' @describeIn LPSOLVE The name of the solver.
#' @export
setMethod("name", "LPSOLVE", function(object) { LPSOLVE_NAME })

#' @describeIn LPSOLVE Imports the Rmosek library.
setMethod("import_solver", "LPSOLVE", function(solver) {
    requireNamespace("lpSolveAPI", quietly = TRUE)
})

setMethod("split_constr", "LPSOLVE", function(solver, constr_map) {
    list(eq_constr = c(constr_map[[EQ_MAP]], constr_map[[LEQ_MAP]]), ineq_constr = list(), nonlin_constr = list())
})

#' @param objective A list representing the canonicalized objective.
#' @param constraints A list of canonicalized constraints.
#' @param cached_data A list mapping solver name to cached problem data.
#' @param warm_start A logical value indicating whether the previous solver result should be used to warm start.
#' @param verbose A logical value indicating whether to print solver output.
#' @param ... Additional arguments to the solver.
#' @describeIn LPSOLVE Call the solver on the canonicalized problem.
setMethod("Solver.solve", "LPSOLVE", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {

    solver_opts <- list(...)
    if (verbose) {
       solver_opts$verbose <- 'normal'
    }

    data <- Solver.get_problem_data(solver, objective, constraints, cached_data)

    c <- data[[C_KEY]]
    dims <- data[[DIMS]]
    nvar <- length(c)
    A <- data[[A_KEY]]
    b <- data[[B_KEY]]
    ncons <- nrow(A)

    lp <- lpSolveAPI::make.lp(ncons, nvar)
    for (i in seq_len(nvar)) {
        lpSolveAPI::set.column(lp, i, A[, i])
    }
    lpSolveAPI::set.rhs(lp, b)
    lpSolveAPI::set.objfn(lp, c)
    lpSolveAPI::set.constr.type(lp, c(rep("=", dims[[EQ_DIM]]),
                                      rep("<=", dims[[LEQ_DIM]]))
                    )
    lpSolveAPI::set.bounds(lp, lower = rep(-Inf, nvar))
    bools <- data[[BOOL_IDX]]
    ints <- data[[INT_IDX]]
    if (length(bools) > 0) {
        lpSolveAPI::set.type(lp, bools, type = "binary")
    }
    if (length(ints) > 0) {
        lpSolveAPI::set.type(lp, ints, type = "integer")
    }

    ## Set the parameters
    do.call(lpSolveAPI::lp.control, c(list(lprec = lp), solver_opts))

    status <- solve(lp)
    results_dict <- list( status, lp)
    names(results_dict) <- c(STATUS, "lprec")

    format_results(solver, results_dict, data, cached_data)
})

#' @param data A list containing information about the problem.
#' @describeIn LPSOLVE Convert raw solver output into standard list of results.
setMethod("format_results", "LPSOLVE", function(solver, results_dict, data, cached_data) {
    new_results <- list()
    new_results[[STATUS]] <- status_map(solver, results_dict[[STATUS]])
    new_results[[SOLVE_TIME]] <- NA
    new_results[[SETUP_TIME]] <- NA
    lp <- results_dict$lprec

    if(new_results[[STATUS]] %in% SOLUTION_PRESENT) {
        ## Get primal variable values
        new_results[[PRIMAL]] <- lpSolveAPI::get.variables(lp)
        new_results[[NUM_ITERS]] <- lpSolveAPI::get.total.iter(lp)
        ## Get objective value
        new_results[[VALUE]] <- lpSolveAPI::get.objective(lp) + data[[OFFSET]]
    }
    new_results
})

