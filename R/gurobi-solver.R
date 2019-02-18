#'
#' The GUROBI class.
#'
#' This class is an interface for the commercial GUROBI solver.
#'
#' @references \emph{Gurobi optimizer reference manual version 5.0,} Gurobi Optimization, Inc., Houston, Texas, July 2012.
#' @seealso the \href{http://www.gurobi.com/documentation/7.5/refman/r_api_overview.html}{GUROBI Official Site}.
#' @name GUROBI-class
#' @aliases GUROBI
#' @rdname GUROBI-class
setClass("GUROBI", contains = "Solver")

#' @rdname GUROBI-class
#' @export
GUROBI <- function() {
  new("GUROBI")
}

#' @param object,solver A \linkS4class{GUROBI} object.
#' @describeIn GUROBI GUROBI can handle linear programs.
setMethod("lp_capable", "GUROBI", function(solver) { TRUE })

#' @describeIn GUROBI GUROBI can handle second-order cone programs.
setMethod("socp_capable", "GUROBI", function(solver) { TRUE })

#' @describeIn GUROBI GUROBI cannot handle semidefinite programs.
setMethod("sdp_capable", "GUROBI", function(solver) { FALSE })

#' @describeIn GUROBI GUROBI cannot handle exponential cone programs.
setMethod("exp_capable", "GUROBI", function(solver) { FALSE })

#' @describeIn GUROBI GUROBI can handle mixed-integer programs.
setMethod("mip_capable", "GUROBI", function(solver) { TRUE })

#'
#' GUROBI Status Map
#'
#' Map of GUROBI status to CVXR status.
#'
#' @param solver A \linkS4class{GUROBI} object.
#' @param status An exit code returned by GUROBI. See the \href{http://www.gurobi.com/documentation/7.5/refman/optimization_status_codes.html}{GUROBI documentation} for details.
#' @return A string indicating the status, either "optimal", "infeasible", "unbounded", "optimal_inaccurate", "infeasible_inaccurate", "unbounded_inaccurate", or "solver_error".
#' @docType methods
#' @rdname GUROBI-status_map
setMethod("status_map", "GUROBI", function(solver, status) {
  if(status == 2)
    OPTIMAL
  else if(status == 3)
    INFEASIBLE
  else if(status == 5)
    UNBOUNDED
  else if(status == 9)
    OPTIMAL_INACCURATE
  else if(status %in% c(4, 6, 7, 8, 10, 11, 12, 13))
    SOLVER_ERROR
  else
    stop("GUROBI status unrecognized: ", status)
})

#' @describeIn GUROBI The name of the solver.
#' @export
setMethod("name", "GUROBI", function(object) { GUROBI_NAME })

#' @describeIn GUROBI Imports the reticulate library to use the python solver.
setMethod("import_solver", "GUROBI", function(solver) {
    if (requireNamespace("reticulate", quietly = TRUE)) {
        cvxr_options <- .CVXR.options
        sp <- np <- gurobiglue <- NULL

        if (reticulate::py_module_available(module = "numpy")) {
            np <- reticulate::import(module = "numpy", delay_load = TRUE)
        }

        if (reticulate::py_module_available(module = "scipy.sparse")) {
            sp <- reticulate::import(module = "scipy.sparse", delay_load = TRUE)
        }

        glue_module_path = system.file("python", package = "CVXR")
        if (reticulate::py_module_available(module = "gurobipy")) {
            ## delay_load removed in reticulate 1.7+
            gurobiglue <- reticulate::import_from_path(module = "gurobiglue",
                                                       glue_module_path)
        }
        cvxr_options$np <- np
        cvxr_options$sp <- sp
        cvxr_options$gurobiglue <- gurobiglue
        assignInMyNamespace(".CVXR.options", cvxr_options)
        !(is.null(np) || is.null(sp) || is.null(gurobiglue))
    } else {
        FALSE
    }
})

setMethod("split_constr", "GUROBI", function(solver, constr_map) {
    list(eq_constr = c(constr_map[[EQ_MAP]], constr_map[[LEQ_MAP]]),
         ineq_constr = list(),
         nonlin_constr = list())
})


#' @param objective A list representing the canonicalized objective.
#' @param constraints A list of canonicalized constraints.
#' @param cached_data A list mapping solver name to cached problem data.
#' @param warm_start A logical value indicating whether the previous solver result should be used to warm start.
#' @param verbose A logical value indicating whether to print solver output.
#' @param ... Additional arguments to the solver.
#' @describeIn GUROBI Call the solver on the canonicalized problem.
setMethod("Solver.solve", "GUROBI", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {

    solver_opts <- list(...)

    data <- Solver.get_problem_data(solver, objective, constraints, cached_data)

    A <- data[[A_KEY]]
    b <- data[[B_KEY]]
    if (length(b) == 1L) b <- list(b)
    ##G <- data[[G_KEY]]
    ##h <- unlist(data[[H_KEY]])
    c <- matrix(unlist(data[[C_KEY]]), ncol = 1)
    offset <- data[[OFFSET]]
    bool_idx <- as.integer(data[[BOOL_IDX]]) - 1L ## Zero-based indices in python
    if (length(bool_idx) == 1L) bool_idx <- list(bool_idx)
    int_idx <- as.integer(data[[INT_IDX]]) -1L ## Zero-based indices in python
    if (length(int_idx) == 1L) int_idx <- list(int_idx)

  dims <- data[[DIMS]]
  if (is.null(dims$s)) {
      dims$s <- list()
  } else {
      dims$s <- as.list(as.integer(dims$s))
  }
  dims$q <- as.list(as.integer(dims$q))
  dims$l <- as.integer(dims$l)
  dims$f <- as.integer(dims$f)
  dims$ep <- as.list(as.integer(dims$ep))

  ## Fix for 0 dimension matrices
    nDims <- dim(A)
    if (prod(nDims) == 0) A <- matrix(nrow = nDims[1], ncol = nDims[2])
  results_dict <- get_gurobiglue()$gurobi_intf(reticulate::r_to_py(A),
                                               b,
                                               c,
                                               bool_idx,
                                               int_idx,
                                               dims,
                                               offset,
                                               reticulate::dict(solver_opts),
                                               verbose)
    results_dict
    ##format_results(solver, results_dict, data, cached_data)
})

