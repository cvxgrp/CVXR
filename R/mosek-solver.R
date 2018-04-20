
#'
#' The MOSEK class.
#'
#' This class is an interface for the commercial MOSEK solver.
#'
#' @references E. Andersen and K. Andersen. "The MOSEK Interior Point Optimizer for Linear Programming: an Implementation of the Homogeneous Algorithm." \emph{High Performance Optimization}, vol. 33, pp. 197-232, 2000.
#' @seealso the \href{https://www.mosek.com/products/mosek/}{MOSEK Official Site}.
#' @name MOSEK-class
#' @aliases MOSEK
#' @rdname MOSEK-class
setClass("MOSEK", contains = "Solver")

#' @rdname MOSEK-class
#' @export
MOSEK <- function() {
  new("MOSEK")
}

#' @param object,solver A \linkS4class{MOSEK} object.
#' @describeIn MOSEK MOSEK can handle linear programs.
setMethod("lp_capable", "MOSEK", function(solver) { TRUE })

#' @describeIn MOSEK MOSEK can handle second-order cone programs.
setMethod("socp_capable", "MOSEK", function(solver) { TRUE })

#' @describeIn MOSEK MOSEK can handle semidefinite programs.
setMethod("sdp_capable", "MOSEK", function(solver) { TRUE })

#' @describeIn MOSEK MOSEK cannot handle exponential cone programs.
setMethod("exp_capable", "MOSEK", function(solver) { FALSE })

#' @describeIn MOSEK MOSEK cannot handle mixed-integer programs.
setMethod("mip_capable", "MOSEK", function(solver) { FALSE })

#'
#' MOSEK Status Map
#'
#' Map of MOSEK status to CVXR status.
#'
#' @param solver A \linkS4class{MOSEK} object.
#' @param status An exit code returned by MOSEK. See the \href{http://docs.mosek.com/8.0/dotnetfusion/solution_status.html}{MOSEK documentation} for details.
#' @return A string indicating the status, either "optimal", "infeasible", "unbounded", "optimal_inaccurate", "infeasible_inaccurate", "unbounded_inaccurate", or "solver_error".
#' @docType methods
#' @rdname MOSEK-status_map
setMethod("status_map", "MOSEK", function(solver, status) {
  if(status == "OPTIMAL")
    OPTIMAL
  else if(status == "PRIMAL_INFEASIBLE_CER")
    INFEASIBLE
  else if(status == "DUAL_INFEASIBLE_CER")
    UNBOUNDED
  else if(status == "NEAR_OPTIMAL")
    OPTIMAL_INACCURATE
  else if(status == "NEAR_PRIMAL_INFEASIBLE_CER")
    INFEASIBLE_INACCURATE
  else if(status == "NEAR_DUAL_INFEASIBLE_CER")
    UNBOUNDED_INACCURATE
  else
    SOLVER_ERROR
})
#
#' @describeIn MOSEK The name of the solver.
#' @export
setMethod("name", "MOSEK", function(object) { MOSEK_NAME })

#' @describeIn MOSEK Imports the Rmosek library.
setMethod("import_solver", "MOSEK", function(solver) {
    if (requireNamespace("reticulate", quietly = TRUE)) {
        cvxr_options <- .CVXR.options

        sp <- np <- mosekglue <- NULL

        if (reticulate::py_module_available(module = "numpy")) {
            np <- reticulate::import(module = "numpy", delay_load = TRUE)
        }

        if (reticulate::py_module_available(module = "scipy.sparse")) {
            sp <- reticulate::import(module = "scipy.sparse", delay_load = TRUE)
        }

        glue_module_path = system.file("python", package="CVXR")
        if (reticulate::py_module_available(module = "mosek")) {
            mosekglue <- reticulate::import_from_path(module = "mosekglue",
                                                      glue_module_path,
                                                      delay_load = TRUE)
        }
        cvxr_options$np <- np
        cvxr_options$sp <- sp
        cvxr_options$mosekglue <- mosekglue
        assignInMyNamespace(".CVXR.options", cvxr_options)
        !(is.null(np) || is.null(sp) || is.null(mosekglue))
    } else {
        FALSE
    }
})

setMethod("split_constr", "MOSEK", function(solver, constr_map) {
  list(eq_constr = constr_map[[EQ_MAP]], ineq_constr = constr_map[[LEQ_MAP]], nonlin_constr = list())
})

#' @param objective A list representing the canonicalized objective.
#' @param constraints A list of canonicalized constraints.
#' @param cached_data A list mapping solver name to cached problem data.
#' @param warm_start A logical value indicating whether the previous solver result should be used to warm start.
#' @param verbose A logical value indicating whether to print solver output.
#' @param ... Additional arguments to the solver.
#' @describeIn MOSEK Call the solver on the canonicalized problem.
setMethod("Solver.solve", "MOSEK", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {

  solver_opts <- list(...)

  data <- Solver.get_problem_data(solver, objective, constraints, cached_data)

  A <- data[[A_KEY]]
  b <- data[[B_KEY]]
  G <- data[[G_KEY]]
  h <- unlist(data[[H_KEY]])
  c <- unlist(data[[C_KEY]])
  offset <- data[[OFFSET]]

  dims <- data[[DIMS]]
  if (is.null(dims$s)) {
      dims$s <- list()
  } else {
      dims$s <- as.integer(dims$s)
  }
  dims$q <- as.integer(dims$q)
  dims$l <- as.integer(dims$l)
  dims$f <- as.integer(dims$f)
  dims$ep <- as.integer(dims$ep)

  results_dict <- get_mosekglue()$mosek_intf(r2py_sparse(A),
                                             b,
                                             r2py_sparse(G),
                                             h,
                                             c,
                                             dims,
                                             offset,
                                             reticulate::dict(solver_opts),
                                             verbose)
    results_dict
    if (is.null(results_dict)) {

    }
    results_dict
    ##format_results(solver, results_dict, data, cached_data)
})

