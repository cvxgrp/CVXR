#' An interface for the DIFFCP solver, a differentiable wrapper of SCS and ECOS
#'
#' @name DIFFCP-class
#' @aliases DIFFCP
#' @rdname DIFFCP-class
#' @export
setClass("DIFFCP", contains = "SCS")

setMethod("initialize", "SCS",
          function(.Object, ...) {
            ##.Object <- callNextMethod(.Object, ...)
            ## Ensure EXP_CONE_ORDER is set
            .Object@EXP_CONE_ORDER <- c(0L, 1L, 2L)
            .Object@SUPPORTED_CONSTRAINTS <- c(.Object@SUPPORTED_CONSTRAINTS, "SOC", "ExpCone", "PSDConstraint")
            .Object@REQUIRES_CONSTRAINTS <- TRUE
            .Object
          })

## Not needed usually
## #' @rdname SCS-class
## #' @export
## SCS <- function() { new("SCS") }

# Solver capabilities.
#' @describeIn SCS Can the solver handle mixed-integer programs?
setMethod("mip_capable", "SCS", function(solver) { FALSE })
##setMethod("requires_constr", "SCS", function(solver) { TRUE })
##setMethod("supported_constraints", "SCS", function(solver) { c(supported_constraints(ConicSolver()), "SOC", "ExpCone", "PSDConstraint") })

# Map of SCS status to CVXR status.
#' @param solver,object,x A \linkS4class{SCS} object.
#' @param status A status code returned by the solver.
#' @describeIn SCS Converts status returned by SCS solver to its respective CVXPY status.
setMethod("status_map", "SCS", function(solver, status) {
  if(status == 1)
    return(OPTIMAL)
  else if(status == 2)
    return(OPTIMAL_INACCURATE)
  else if(status == -1)
    return(UNBOUNDED)
  else if(status == -6)
    return(UNBOUNDED_INACCURATE)
  else if(status == -2)
    return(INFEASIBLE)
  else if(status == -7)
    return(INFEASIBLE_INACCURATE)
  else if(status %in% c(-3, -4, -5))
    return(SOLVER_ERROR)
  else
    stop("DIFFCP status unrecognized: ", status)
})

#' @describeIn DIFFCP Returns the name of the solver
setMethod("name", "DIFFCP", function(x) { DIFFCP_NAME })

#' @describeIn DIFFCP Imports the solver
##setMethod("import_solver", "DIFFCP", function(solver) { requireNamespace("scs", quietly = TRUE) })
setMethod("import_solver", "DIFFCP", function(solver) { TRUE }) ## we require scs

#' @describeIn DIFFCP supports quadratic objective?
setMethod("supports_quad_obj", "DIFFCP", function(solver) { FALSE }) ## This one doesn't

#' @describeIn DIFFCP Returns a new problem and data for inverting the new solution
setMethod("perform", signature(object = "DIFFCP", problem = "Problem"), function(object, problem) {
  results <- ConicSolver.prepare_data_and_inv_data(problem)
  mats <- apply_parameters(results$problem, keep_zeros = TRUE)
  # Returns a new problem and data for inverting the new solution.
  data <- results$data
  data[[C_KEY]] <- mats$c
  inv_data <- results$inv_data
  inv_data[[OFFSET_KEY]] <- mats$d
  data[[A_KEY]] <- mats$A
  data[[B_KEY]] <- mats$b
  list(data = data, inv_data = inv_data)
}

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A list containing data necessary for the inversion.
#' @describeIn DIFFCP Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "DIFFCP", solution = "list", inverse_data = "list"), function(object, solution, inverse_data) {

  # Returns the solution to the original problem given the inverse_data.
  status <- status_map(object, solution$info$status)

  attr <- list()
  attr[[SOLVE_TIME]] <- solution$info$solveTime
  attr[[SETUP_TIME]] <- solution$info$setupTime
  attr[[NUM_ITERS]] <- solution$info$iter
  attr[[EXTRA_STATS]] <- solution

  if(status %in% SOLUTION_PRESENT) {
    primal_val <- solution$info$pobj
    opt_val <- primal_val + inverse_data[[OFFSET]]
    primal_vars <- list()
    var_id <- inverse_data[[object@var_id]]
    primal_vars[[as.character(var_id)]] <- as.matrix(solution$x)

    num_zero <- inverse_data[[object@dims]]@zero
    eq_idx <- seq_len(num_zero)
    ineq_idx <- seq(num_zero + 1, length.out = length(solution$y) - num_zero)
    eq_dual_vars <- get_dual_values(solution$y[eq_idx], DIFFCP.extract_dual_value, inverse_data[[object@eq_constr]])
    ineq_dual_vars <- get_dual_values(solution$y[ineq_idx], DIFFCP.extract_dual_value, inverse_data[[object@neq_constr]])

    dual_vars <- list()
    dual_vars <- utils::modifyList(dual_vars, eq_dual_vars)
    dual_vars <- utils::modifyList(dual_vars, ineq_dual_vars)
    return(Solution(status, opt_val, primal_vars, dual_vars, attr))
  } else
    return(failure_solution(status))
})

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance on the primal and dual residual.
#' @param reltol The relative tolerance on the duality gap.
#' @param abstol The absolute tolerance on the duality gap.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn DIFFCP Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "DIFFCP", function(object, data, warm_start, verbose, feastol, reltol, abstol,
                                            num_iter, solver_opts, solver_cache = new.env(parent = emptyenv()) ) {
  A  <- data[[A_KEY]]
  if (!inherits(A, "dgCMatrix")) A  <- as(as(A, "CsparseMatrix"), "generalMatrix")

  args <- list(A = A, b = data[[B_KEY]], c = data[[C_KEY]])
  if(warm_start && !is.null(solver_cache) && length(solver_cache) > 0 && name(object) %in% names(solver_cache)) {
    obj_name <- name(object)
    args$x <- solver_cache[[obj_name]]$x
    args$y <- solver_cache[[obj_name]]$y
    args$s <- solver_cache[[obj_name]]$s
  }
  cones <- DIFFCP.dims_to_solver_dict(data[[ConicSolver()@dims]])

  ## if(!all(c(is.null(feastol), is.null(reltol), is.null(abstol)))) {
  ##   warning("Ignoring inapplicable parameter feastol/reltol/abstol for DIFFCP.")
  ## }
  solver_defaults  <- SOLVER_DEFAULT_PARAM$DIFFCP

  if(is.null(num_iter)) {
    num_iter <- solver_defaults$max_iters
  }
  if (is.null(reltol)) {
    reltol <- solver_defaults$eps_rel
  }
  if (is.null(abstol)) {
    abstol  <- solver_defaults$eps_abs
  }
  if (is.null(feastol)) {
    feastol  <- solver_defaults$eps_infeas
  }

  control =  scs::scs_control(max_iters = num_iter, verbose = verbose, eps_rel = reltol, eps_abs = abstol, eps_infeas = feastol)

  #Fill in parameter values
  control[names(solver_opts)] <- solver_opts

  # Returns the result of the call to the solver.
  results <- scs::scs(A = args$A, b = args$b, obj = args$c, cone = cones, control = control)
  if(!is.null(solver_cache) && length(solver_cache) > 0)
    solver_cache[[name(object)]] <- results
  return(results)
})
