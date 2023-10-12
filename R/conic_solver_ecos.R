#'
#' An interface for the ECOS solver
#'
#' @name ECOS-class
#' @aliases ECOS
#' @rdname ECOS-class
#' @export
setClass("ECOS", representation(exp_cone_order = "integer"),   # Order of exponential cone arguments for solver. Internal only!
                 prototype(exp_cone_order = c(0L, 2L, 1L)), contains = "ConicSolver")

#' @rdname ECOS-class
#' @export
ECOS <- function() { new("ECOS") }

# Solver capabilities.
#' @describeIn ECOS Can the solver handle mixed-integer programs?
setMethod("mip_capable", "ECOS", function(solver) { FALSE })
setMethod("supported_constraints", "ECOS", function(solver) { c(supported_constraints(ConicSolver()), "SOC", "ExpCone") })

# EXITCODES from ECOS
# ECOS_OPTIMAL  (0)   Problem solved to optimality
# ECOS_PINF     (1)   Found certificate of primal infeasibility
# ECOS_DINF     (2)   Found certificate of dual infeasibility
# ECOS_INACC_OFFSET (10)  Offset exitflag at inaccurate results
# ECOS_MAXIT    (-1)  Maximum number of iterations reached
# ECOS_NUMERICS (-2)  Search direction unreliable
# ECOS_OUTCONE  (-3)  s or z got outside the cone, numerics?
# ECOS_SIGINT   (-4)  solver interrupted by a signal/ctrl-c
# ECOS_FATAL    (-7)  Unknown problem in solver

# Map of ECOS status to CVXR status.
#' @param solver,object,x A \linkS4class{ECOS} object.
#' @param status A status code returned by the solver.
#' @describeIn ECOS Converts status returned by the ECOS solver to its respective CVXPY status.
setMethod("status_map", "ECOS", function(solver, status) {
  if(status == 0)
    return(OPTIMAL)
  else if(status == 1)
    return(INFEASIBLE)
  else if(status == 2)
    return(UNBOUNDED)
  else if(status == 10)
    return(OPTIMAL_INACCURATE)
  else if(status == 11)
    return(INFEASIBLE_INACCURATE)
  else if(status == 12)
    return(UNBOUNDED_INACCURATE)
  else if(status %in% c(-1, -2, -3, -4, -7))
    return(SOLVER_ERROR)
  else
    stop("ECOS status unrecognized: ", status)
})

#' @describeIn ECOS Imports the solver
##setMethod("import_solver", "ECOS", function(solver) { requireNamespace("ECOSolveR", quietly = TRUE) })
## Since ECOS is required, this is always TRUE
setMethod("import_solver", "ECOS", function(solver) { TRUE })

#' @describeIn ECOS Returns the name of the solver
setMethod("name", "ECOS", function(x) { ECOS_NAME })

#' @param problem A \linkS4class{Problem} object.
#' @describeIn ECOS Returns a new problem and data for inverting the new solution.
setMethod("perform", signature(object = "ECOS", problem = "Problem"), function(object, problem) {
  ## SHOULD problem not be of class ParamConeProg ?
  if (! problem@formatted) {
    problem <- format_constr(object, problem, problem@exp_cone_order)
  }
  data <- list(PARAM_PROB = problem)
  inv_data <- list()
  data[[object@dims]] <- inv_data[[object@dims]] <- problem@cone_dims

  constr_map <- problem@constr_map
  inv_data[[EQ_CONSTR]] <- constr_map$ZeroConstraint

  inv_data[[object@dims]] <- id(variables(problem)[[1L]])
  offsets <- ConicSolver.get_coeff_offset(problem@objective@args[[1L]])

  data[[C_KEY]] <- as.vector(offsets[[1]])
  data[[OFFSET]] <- offsets[[2]]
  inv_data[[OFFSET]] <- data[[OFFSET]][1]
  inv_data[[object@eq_constr]] <- constr_map$ZeroConstraint
  offsets <- group_coeff_offset(object, problem, constr_map$ZeroConstraint, ECOS()@exp_cone_order)
  data[[A_KEY]] <- offsets[[1]]
  data[[B_KEY]] <- offsets[[2]]

  # Order and group nonlinear constraints.
  neq_constr <- c(constr_map$NonPosConstraint, constr_map$SOC, constr_map$ExpCone)
  inv_data[[object@neq_constr]] <- neq_constr
  offsets <- group_coeff_offset(object, problem, neq_constr, ECOS()@exp_cone_order)
  data[[G_KEY]] <- offsets[[1]]
  data[[H_KEY]] <- offsets[[2]]

  return(list(object, data, inv_data))
})

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A list containing data necessary for the inversion.
#' @describeIn ECOS Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "ECOS", solution = "list", inverse_data = "list"), function(object, solution, inverse_data) {
  status <- status_map(object, solution$retcodes[["exitFlag"]])

  # Timing data.
  attr <- list()
  attr[[SOLVE_TIME]] <- solution$timing[["tsolve"]]
  attr[[SETUP_TIME]] <- solution$timing[["tsetup"]]
  attr[[NUM_ITERS]] <- solution$retcodes[["iter"]]

  if(status %in% SOLUTION_PRESENT) {
    primal_val <- solution$summary[["pcost"]]
    opt_val <- primal_val + inverse_data[[OFFSET]]
    primal_vars <- list()
    var_id <- inverse_data[[object@var_id]]
    primal_vars[[as.character(var_id)]] <- as.matrix(solution$x)

    eq_dual <- get_dual_values(solution$y, extract_dual_value, inverse_data[[object@eq_constr]])
    leq_dual <- get_dual_values(solution$z, extract_dual_value, inverse_data[[object@neq_constr]])
    eq_dual <- utils::modifyList(eq_dual, leq_dual)
    dual_vars <- eq_dual

    return(Solution(status, opt_val, primal_vars, dual_vars, attr))
  } else
    return(failure_solution(status))
})

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose An integer number indicating level of solver verbosity.
#' @param feastol The feasible tolerance on the primal and dual residual.
#' @param reltol The relative tolerance on the duality gap.
#' @param abstol The absolute tolerance on the duality gap.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn ReductionSolver Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "ECOS", function(object, data, warm_start, verbose,
                                             feastol,
                                             reltol,
                                             abstol,
                                             num_iter,
                                             solver_opts, solver_cache) {

  if (missing(solver_cache)) solver_cache <- new.env(parent=emptyenv())

  cones <- ECOS.dims_to_solver_dict(data[[ConicSolver()@dims]])

  if(is.null(feastol)) {
      feastol <- SOLVER_DEFAULT_PARAM$ECOS$feastol
  }
  if(is.null(reltol)) {
      reltol <- SOLVER_DEFAULT_PARAM$ECOS$reltol
  }
  if(is.null(abstol)) {
      abstol <- SOLVER_DEFAULT_PARAM$ECOS$abstol
  }
  if(is.null(num_iter)) {
      num_iter <- SOLVER_DEFAULT_PARAM$ECOS$maxit
  }

  ecos_opts <- ECOSolveR::ecos.control(maxit = as.integer(num_iter), feastol = feastol, reltol = reltol, abstol = abstol, verbose = as.integer(verbose))
  ecos_opts[names(solver_opts)] <- solver_opts
  solution <- ECOSolveR::ECOS_csolve(c = data[[C_KEY]], G = data[[G_KEY]], h = data[[H_KEY]], dims = cones, A = data[[A_KEY]], b = data[[B_KEY]], control = ecos_opts)
  return(solution)
})
