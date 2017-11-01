#'
#' The Minimize class.
#'
#' This class represents an optimization objective for minimization.
#'
#' @slot expr A scalar \linkS4class{Expression} to minimize.
#' @rdname Minimize
#' @export
Minimize <- setClass("Minimize", representation(expr = "ConstValORExpr"), contains = "Canonical")

setMethod("initialize", "Minimize", function(.Object, expr) {
    .Object@expr <- as.Constant(expr)
    if(!all(size(.Object@expr) == c(1,1)))
      stop("The objective must resolve to a scalar")
    return(.Object)
})

#' @rdname canonicalize
#' @describeIn Minimize Pass on the target expression's objective and constraints.
setMethod("canonicalize", "Minimize", function(object) { canonical_form(object@expr) })

#' @rdname expression-parts
#' @describeIn Minimize Returns the \linkS4class{Variable} objects in the objective.
setMethod("variables", "Minimize", function(object) { variables(object@expr) })

#' @rdname expression-parts
#' @describeIn Minimize Returns the \linkS4class{Parameter} objects in the objective.
setMethod("parameters", "Minimize", function(object) { parameters(object@expr) })

#' @rdname expression-parts
#' @describeIn Minimize Returns the \linkS4class{Constant} objects in the objective.
setMethod("constants", "Minimize", function(object) { constants(object@expr) })

#' @rdname is_dcp
#' @describeIn Minimize A logical value indicating whether the objective is convex.
setMethod("is_dcp", "Minimize", function(object) { is_convex(object@expr) })

#' @rdname value-methods
#' @describeIn Minimize The value of the objective expression.
setMethod("value", "Minimize", function(object) { value(object@expr) })

# The value of the objective given the solver primal value.
setMethod("primal_to_result", "Minimize", function(object, result) { result })

#'
#' The Maximize class.
#'
#' This class represents an optimization objective for maximization.
#'
#' @slot expr A scalar \linkS4class{Expression} to maximize.
#' @examples
#' x <- Variable(3)
#' alpha <- c(0.8,1.0,1.2)
#' obj <- sum(log(alpha + x))
#' constr <- list(x >= 0, sum(x) == 1)
#' prob <- Problem(Maximize(obj), constr)
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' @rdname Maximize
#' @export
Maximize <- setClass("Maximize", contains = "Minimize")

#'
#' Arithmetic Operations on Objectives
#' 
#' Add, subtract, multiply, or divide optimization objectives.
#' 
#' @param e1 The left-hand \linkS4class{Minimize}, \linkS4class{Maximize}, or numeric value.
#' @param e2 The right-hand \linkS4class{Minimize}, \linkS4class{Maximize}, or numeric value.
#' @return A \linkS4class{Minimize} or \linkS4class{Maximize} object.
#' @name Objective-arith
NULL

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Minimize", e2 = "missing"), function(e1, e2) { Maximize(expr = -e1@expr) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Minimize", e2 = "Minimize"), function(e1, e2) { Minimize(e1@expr + e2@expr) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Minimize", e2 = "Maximize"), function(e1, e2) { stop("Problem does not follow DCP rules") })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Minimize", e2 = "numeric"), function(e1, e2) { if(length(e2) == 1 && e2 == 0) e1 else stop("Unimplemented") })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "numeric", e2 = "Minimize"), function(e1, e2) { e2 + e1 })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Minimize", e2 = "Minimize"), function(e1, e2) { e1 + (-e2) })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Minimize", e2 = "Maximize"), function(e1, e2) { e1 + (-e2) })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Minimize", e2 = "numeric"), function(e1, e2) { e1 + (-e2) })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "numeric", e2 = "Minimize"), function(e1, e2) { if(length(e2) == 1 && e2 == 0) -e1 else stop("Unimplemented") })

#' @rdname Objective-arith
setMethod("*", signature(e1 = "Minimize", e2 = "numeric"), function(e1, e2) {
  if(e2 >= 0) Minimize(expr = e1@expr * e2) else Maximize(expr = e1@expr * e2)
})

#' @rdname Objective-arith
setMethod("*", signature(e1 = "Maximize", e2 = "numeric"), function(e1, e2) {
  if(e2 < 0) Minimize(expr = e1@expr * e2) else Maximize(expr = e1@expr * e2)
})

#' @rdname Objective-arith
setMethod("*", signature(e1 = "numeric", e2 = "Minimize"), function(e1, e2) { e2 * e1 })

#' @rdname Objective-arith
setMethod("/", signature(e1 = "Minimize", e2 = "numeric"), function(e1, e2) { e1 * (1.0/e2) })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Maximize", e2 = "missing"), function(e1, e2) { Minimize(expr = -e1@expr) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Maximize", e2 = "Maximize"), function(e1, e2) { Maximize(expr = e1@expr + e2@expr) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Maximize", e2 = "Minimize"), function(e1, e2) { stop("Problem does not follow DCP rules") })

#' @rdname canonicalize
#' @describeIn Maximize Negates the target expression's objective.
setMethod("canonicalize", "Maximize", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constraints <- canon[[2]]
  list(lo.neg_expr(obj), constraints)
})

#' @rdname is_dcp
#' @describeIn Maximize A logical value indicating whether the objective is concave.
setMethod("is_dcp", "Maximize", function(object) { is_concave(object@expr) })

#' @rdname curvature-methods
#' @describeIn Maximize A logical value indicating whether the objective is quadratic.
setMethod("is_quadratic", "Maximize", function(object) { is_quadratic(object@expr) })

# The value of the objective given the solver primal value.
setMethod("primal_to_result", "Maximize", function(object, result) { -result })

#'
#' The SolverStats class.
#' 
#' This class contains the miscellaneous information that is returned by a solver after solving, but that is not captured directly by the \linkS4class{Problem} object.
#' 
#' @slot solver_name The name of the solver.
#' @slot solve_time The time (in seconds) it took for the solver to solve the problem.
#' @slot setup_time The time (in seconds) it took for the solver to set up the problem.
#' @slot num_iters The number of iterations the solver had to go through to find a solution.
#' @name SolverStats-class
#' @rdname SolverStats-class
#' @export
.SolverStats <- setClass("SolverStats", representation(solver_name = "character", solve_time = "numeric", setup_time = "numeric", num_iters = "numeric"),
                         prototype(solver_name = NA_character_, solve_time = NA_real_, setup_time = NA_real_, num_iters = NA_real_))

#'
#' Solver Statistics
#'
#' Reports some of the miscellaneous information that is returned by a solver after solving, but that is not captured directly by the \linkS4class{Problem} object.
#' 
#' @param results_dict A list containing the results returned by the solver.
#' @param solver_name The name of the solver.
#' @return A list containing
#' \itemize{
#'   \item{solver_name}{The name of the solver.}
#'   \item{solve_time}{The time (in seconds) it took for the solver to solve the problem.}
#'   \item{setup_time}{The time (in seconds) it took for the solver to set up the problem.}
#'   \item{num_iters}{The number of iterations the solver had to go through to find a solution.}
#' }
#' @name SolverStats
#' @rdname SolverStats-class
#' @export
SolverStats <- function(results_dict = list(), solver_name = NA_character_) {
    solve_time <- NA_real_
    setup_time <- NA_real_
    num_iters <- NA_real_

    if(SOLVE_TIME %in% names(results_dict))
        solve_time <- results_dict[[SOLVE_TIME]]
    if(SETUP_TIME %in% names(results_dict))
        setup_time <- results_dict[[SETUP_TIME]]
    if(NUM_ITERS %in% names(results_dict))
        num_iters <- results_dict[[NUM_ITERS]]
    list(SOLVER = solver_name,
         SOLVE_TIME = solve_time,
         SETUP_TIME = setup_time,
         NUM_ITERS = num_iters)
    ##  .SolverStats(solver_name = solver_name, solve_time = solve_time, setup_time = setup_time, num_iters = num_iters)
}

#'
#' The SizeMetrics class.
#'
#' This class contains various metrics regarding the problem size.
#' 
#' @slot num_scalar_variables The number of scalar variables in the problem.
#' @slot num_scalar_data The number of constants used across all matrices and vectors in the problem. Some constants are not apparent when the problem is constructed. For example, the \code{sum_squares} expression is a wrapper for a \code{quad_over_lin} expression with a constant \code{1} in the denominator.
#' @slot num_scalar_eq_constr The number of scalar equality constraints in the problem.
#' @slot num_scalar_leq_constr The number of scalar inequality constraints in the problem.
#' @slot max_data_dimension The longest dimension of any data block constraint or parameter.
#' @slot max_big_small_squared The maximum value of (big)(small)^2 over all data blocks of the problem, where (big) is the larger dimension and (small) is the smaller dimension for each data block.
#' @name SizeMetrics-class
#' @rdname SizeMetrics-class
#' @export
.SizeMetrics <- setClass("SizeMetrics", representation(num_scalar_variables = "numeric", num_scalar_data = "numeric", num_scalar_eq_constr = "numeric", num_scalar_leq_constr = "numeric",
                                                       max_data_dimension = "numeric", max_big_small_squared = "numeric"),
                         prototype(num_scalar_variables = NA_real_, num_scalar_data = NA_real_, num_scalar_eq_constr = NA_real_, num_scalar_leq_constr = NA_real_,
                                   max_data_dimension = NA_real_, max_big_small_squared = NA_real_))

#'
#' Problem Size Metrics
#'
#' Reports various metrics regarding the problem size.
#' 
#' @param problem A \linkS4class{Problem} object.
#' @return A \linkS4class{SizeMetrics} object containing size metrics of the input problem.
#' @name SizeMetrics
#' @rdname SizeMetrics-class
#' @export
SizeMetrics <- function(problem) {
  # num_scalar_variables
  num_scalar_variables <- 0
  for(var in variables(problem))
    num_scalar_variables <- num_scalar_variables + prod(size(var))

  # num_scalar_data, max_data_dimension, and max_big_small_squared
  max_data_dimension <- 0
  num_scalar_data <- 0
  max_big_small_squared <- 0
  for(const in c(constants(problem), parameters(problem))) {
    big <- 0
    # Compute number of data
    num_scalar_data <- num_scalar_data + prod(size(const))
    big <- max(size(const))
    small <- min(size(const))

    # Get max data dimension
    if(max_data_dimension < big)
      max_data_dimension <- big

    if(max_big_small_squared < as.numeric(big)*small*small)
      max_big_small_squared <- as.numeric(big)*small*small
  }

  # num_scalar_eq_constr
  num_scalar_eq_constr <- 0
  for(constraint in problem@constraints) {
    if(is(constraint, "EqConstraint"))
      num_scalar_eq_constr <- num_scalar_eq_constr + prod(size(constraint@.expr))
  }

  # num_scalar_leq_constr
  num_scalar_leq_constr <- 0
  for(constraint in problem@constraints) {
    if(is(constraint, "LeqConstraint") && class(constraint) == "LeqConstraint")
      num_scalar_leq_constr <- num_scalar_leq_constr + prod(size(constraint@.expr))
  }

  .SizeMetrics(num_scalar_variables = num_scalar_variables, num_scalar_data = num_scalar_data, num_scalar_eq_constr = num_scalar_eq_constr, num_scalar_leq_constr = num_scalar_leq_constr,
               max_data_dimension = max_data_dimension, max_big_small_squared = max_big_small_squared)
}

setClassUnion("SizeMetricsORNull", c("SizeMetrics", "NULL"))

#'
#' The Problem class.
#'
#' This class represents a convex optimization problem.
#'
#' @slot objective A \linkS4class{Minimize} or \linkS4class{Maximize} object representing the optimization objective.
#' @slot constraints (Optional) A list of constraints on the optimization variables.
#' @slot value (Internal) Used internally to hold the value of the optimization objective at the solution.
#' @slot status (Internal) Used internally to hold the status of the problem solution.
#' @slot .cached_data (Internal) Used internally to hold cached matrix data.
#' @slot .separable_problems (Internal) Used internally to hold separable problem data.
#' @slot .size_metrics (Internal) Used internally to hold size metrics.
#' @slot .solver_stats (Internal) Used internally to hold solver statistics.
#' @name Problem-class
#' @rdname Problem-class
#' @export
.Problem <- setClass("Problem", representation(objective = "Minimize", constraints = "list", value = "numeric", status = "character", .cached_data = "list", .separable_problems = "list", .size_metrics = "SizeMetricsORNull", .solver_stats = "list"),
                    prototype(constraints = list(), value = NA_real_, status = NA_character_, .cached_data = list(), .separable_problems = list(), .size_metrics = NULL, .solver_stats = NULL),
                    validity = function(object) {
                      if(!(class(object@objective) %in% c("Minimize", "Maximize")))
                        stop("[Problem: objective] objective must be of class Minimize or Maximize")
                      if(!is.na(object@value))
                        stop("[Problem: value] value should not be set by user")
                      if(!is.na(object@status))
                        stop("[Problem: status] status should not be set by user")
                      if(length(object@.cached_data) > 0)
                        stop("[Problem: .cached_data] .cached_data is an internal slot and should not be set by user")
                      if(length(object@.separable_problems) > 0)
                        stop("[Problem: .separable_problems] .separable_problems is an internal slot and should not be set by user")
                      if(length(object@.solver_stats) > 0)
                        stop("[Problem: .solver_stats] .solver_stats is an internal slot and should not be set by user")
                      return(TRUE)
                    }, contains = "Canonical")

#'
#' Problem Constructor
#'
#' Construct a \linkS4class{Problem} object.
#' 
#' @param objective A \linkS4class{Minimize} or \linkS4class{Maximize} object representing the optimization objective.
#' @param constraints (Optional) A list of constraints on the optimization variables.
#' @return A \linkS4class{Problem} object.
#' @name Problem
#' @rdname Problem-class
#' @export
Problem <- function(objective, constraints = list()) {
  .Problem(objective = objective, constraints = constraints)
}

# Used in problem@.cached_data to check if the problem's objective or constraints have changed.
CachedProblem <- function(objective, constraints) { list(objective = objective, constraints = constraints, class = "CachedProblem") }

# Used by pool.map to send solve result back. Unsure if this is necessary for multithreaded operation in R.
SolveResult <- function(opt_value, status, primal_values, dual_values) { list(opt_value = opt_value, status = status, primal_values = primal_values, dual_values = dual_values, class = "SolveResult") }

#'
#' Problem Initialization
#'
#' @name Problem
#' @rdname Problem-class
setMethod("initialize", "Problem", function(.Object, ..., objective, constraints = list(), value = NA_real_, status = NA_character_, .cached_data = list(), .separable_problems = list(), .size_metrics = SizeMetrics(), .solver_stats = list()) {
  .Object@objective <- objective
  .Object@constraints <- constraints
  .Object@value <- value
  .Object@status <- status

  # Cached processed data for each solver.
  .reset_cache <- function(object) {
    for(solver_name in names(SOLVERS))
      object@.cached_data[[solver_name]] <- ProblemData()
    object@.cached_data[[PARALLEL]] <- CachedProblem(NA, NULL)
    object
  }
  .Object@.cached_data <- .cached_data
  .Object <- .reset_cache(.Object)

  # List of separable (sub)problems
  .Object@.separable_problems <- .separable_problems

  # Information about the size of the problem and its constituent parts
  .Object@.size_metrics <- SizeMetrics(.Object)

  # Benchmarks reported by the solver
  .Object@.solver_stats <- .solver_stats
  .Object
})

#' @rdname problem-parts
#' @describeIn Problem The objective of the problem.
setMethod("objective", "Problem", function(object) { object@objective })

#' @rdname problem-parts
#' @describeIn Problem A list of the constraints of the problem.
setMethod("constraints", "Problem", function(object) { object@constraints })

#' @rdname value-methods
#' @describeIn Problem The value from the last time the problem was solved.
setMethod("value", "Problem", function(object) { object@value })

#' @rdname value-methods
#' @describeIn Problem Set the value of optimal objective.
setMethod("value<-", "Problem", function(object, value ) {
    object@value <- value
    object
})

#' @rdname problem-parts
#' @describeIn Problem The status from the last time the problem was solved.
setMethod("status", "Problem", function(object) { object@status })

# Set the status of the problem.
setMethod("status<-", "Problem", function(object, value ) {
    object@status <- value
    object
})

#' @rdname is_dcp
#' @describeIn Problem A logical value indicating whether the problem statisfies DCP rules.
setMethod("is_dcp", "Problem", function(object) {
  all(sapply(c(object@constraints, list(object@objective)), is_dcp))
})

#' @rdname is_qp
#' @describeIn Problem A logical value indicating whether the problem is a quadratic program.
setMethod("is_qp", "Problem", function(object) {
  for(c in object@constraints) {
    if(!(is(c, "EqConstraint") || is_pwl(c@expr)))
      return(FALSE)
  }
  return(is_dcp(object) && is_quadratic(object@objective@expr))
})

#' @rdname canonicalize
#' @describeIn Problem The graph implementation of the problem.
setMethod("canonicalize", "Problem", function(object) {
  obj_canon <- canonical_form(object@objective)
  canon_constr <- obj_canon[[2]]

  for(constr in object@constraints)
    canon_constr <- c(canon_constr, canonical_form(constr)[[2]])
  list(obj_canon[[1]], canon_constr)
})

#' @rdname expression-parts
#' @describeIn Problem List of \linkS4class{Variable} objects in the problem.
setMethod("variables", "Problem", function(object) {
  vars_ <- variables(object@objective)
  constrs_ <- lapply(object@constraints, function(constr) { variables(constr) })
  unique(flatten_list(c(vars_, constrs_)))   # Remove duplicates
})

#' @rdname expression-parts
#' @describeIn Problem List of \linkS4class{Parameter} objects in the problem.
setMethod("parameters", "Problem", function(object) {
  params <- parameters(object@objective)
  constrs_ <- lapply(object@constraints, function(constr) { parameters(constr) })
  unique(flatten_list(c(params, constrs_)))   # Remove duplicates
})

#' @rdname expression-parts
#' @describeIn Problem List of \linkS4class{Constant} objects in the problem.
setMethod("constants", "Problem", function(object) {
  constants_ <- lapply(object@constraints, function(constr) { constants(constr) })
  constants_ <- c(constants(object@objective), constants_)
  unique(flatten_list(constants_))   # TODO: Check duplicated constants are removed correctly
})

#' @rdname problem-parts
#' @describeIn Problem Information about the size of the problem.
setMethod("size_metrics", "Problem", function(object) { object@.size_metrics })

#' @rdname problem-parts
#' @describeIn Problem Additional information returned by the solver.
setMethod("solver_stats", "Problem", function(object) { object@.solver_stats })

# Set the additional information returned by the solver in the problem.
setMethod("solver_stats<-", "Problem", function(object, value ) {
    object@.solver_stats <- value
    object
})

# solve.Problem <- function(object, method, ...) {
#  if(missing(method))
#    .solve(object, ...)
#  else {
#    func <- Problem.REGISTERED_SOLVE_METHODS[func_name]
#    func(object, ...)
#  }
# }

# Problem.register_solve <- function(name, func) {
#  Problem.REGISTERED_SOLVE_METHODS[name] <- func
# }

#' @rdname get_problem_data
setMethod("get_problem_data", signature(object = "Problem", solver = "character"), function(object, solver) {
  canon <- canonicalize(object)
  objective <- canon[[1]]
  constraints <- canon[[2]]

  # Raise an error if the solver cannot handle the problem
  validate_solver(SOLVERS[[solver]], constraints)
  Solver.get_problem_data(SOLVERS[[solver]], objective, constraints, object@.cached_data)
})

#'
#' Solve a DCP Problem
#' 
#' Solve a DCP compliant optimization problem.
#' 
#' @param object A \linkS4class{Problem} object.
#' @param solver (Optional) A string indicating the solver to use. Defaults to "ECOS".
#' @param ignore_dcp (Optional) A logical value indicating whether to override the DCP check for a problem.
#' @param warm_start (Optional) A logical value indicating whether the previous solver result should be used to warm start.
#' @param verbose (Optional) A logical value indicating whether to print additional solver output.
#' @param parallel (Optional) A logical value indicating whether to solve in parallel if the problem is separable.
#' @param ... Additional options that will be passed to the specific solver. In general, these options will override any default settings imposed by CVXR.
#' @return A list containing the optimal value, primal, and dual variables for the problem.
#' @name solve
#' @rdname solve
#' @export
solve.Problem <- function(object, solver, ignore_dcp = FALSE, warm_start = FALSE, verbose = FALSE, parallel = FALSE, ...) {
  if(!is_dcp(object)) {
    if(ignore_dcp)
      print("Problem does not follow DCP rules. Solving a convex relaxation.")
    else
      stop("Problem does not follow DCP rules.")
  }

  # Standard cone problem
  canon <- canonicalize(object)
  objective <- canon[[1]]
  constraints <- canon[[2]]

  # Solve in parallel
  if(parallel) {
    # Check if the objective or constraint has changed
    if(objective != object@.cached_data[[PARALLEL]]@objective ||
       constraints != object@.cached_data[[PARALLEL]]@constraints) {
      object@.separable_problems <- get_separable_problems(object)
      object@.cached_data[[PARALLEL]] <- CachedProblem(objective, constraints)
    }

    if(length(object@.separable_problems) > 1)
      .parallel_solve(object, solver, ignore_dcp, warm_start, verbose, ...)
  }

  # Choose a solver/check the chosen solver
  if(missing(solver) || is.null(solver) || is.na(solver))
    solver <- Solver.choose_solver(constraints)
  else if(solver %in% names(SOLVERS)) {
    solver <- SOLVERS[[solver]]
    validate_solver(solver, constraints)
  } else
    stop("Unknown solver")

  cached_data <- get_sym_data(solver, objective, constraints, object@.cached_data)
  object@.cached_data <- cached_data
  sym_data <- cached_data[[name(solver)]]@sym_data

  # Presolve couldn't solve the problem.
  if(is.na(sym_data@.presolve_status))
    results_dict <- Solver.solve(solver, objective, constraints, object@.cached_data, warm_start, verbose, ...)
  # Presolve determined problem was unbounded or infeasible
  else {
    results_dict <- list(sym_data@.presolve_status)
    names(results_dict) <- STATUS
  }

  # Correct optimal value if the objective is Maximize
  if(results_dict[[STATUS]] %in% SOLUTION_PRESENT)
    results_dict[[VALUE]] <- primal_to_result(object@objective, results_dict[[VALUE]])
  ##return(results_dict)

  ##Update object
  ##status(object) <- results_dict[[STATUS]]
  ##solver_stats(object) <- SolverStats(results_dict, name(solver))

  ## TODO: For now, don't update and save problem state since R can't modify in place variable slots.
  ##object <- .update_problem_state(object, results_dict, sym_data, solver)
  ##return(object)
  return(valuesById(object, results_dict, sym_data, solver))
}

# TODO: Finish implementation of parallel solve.
.parallel_solve.Problem <- function(object, solver = NULL, ignore_dcp = FALSE, warm_start = FALSE, verbose = FALSE, ...) {
  .solve_problem <- function(problem) {
    result <- solve(problem, solver = solver, ignore_dcp = ignore_dcp, warm_start = warm_start, verbose = verbose, parallel = FALSE, ...)
    opt_value <- result$optimal_value
    status <- result$status
    primal_values <- result$primal_values
    dual_values <- result$dual_values
    return(SolveResults(opt_value, status, primal_values, dual_values))
  }

  # TODO: Finish implementation of parallel solve
  statuses <- sapply(solve_results, function(solve_result) { solve_result$status })

  # Check if at least one subproblem is infeasible or inaccurate
  for(status in INF_OR_UNB) {
    if(status %in% statuses) {
      .handle_no_solution(object, status)
      break
    } else {
      for(i in 1:length(solve_results)) {
        subproblem <- object@separable_problems[i]
        solve_result <- solve_results[i]

        for(j in 1:length(solve_result$primal_values)) {
          var <- variables(subproblem)[j]
          primal_value <- solve_result$primal_values[j]
          subproblem <- save_value(var, primal_value)   # TODO: Fix this since R makes copies
        }

        for(j in 1:length(solve_result$dual_values)) {
          constr <- subproblem@constraints[j]
          dual_value <- solve_result$dual_values[j]
          subproblem <- save_value(constr, dual_value)   # TODO: Fix this since R makes copies
        }
      }

      object@value <- sum(sapply(solve_results, function(solve_result) { solve_result$optimal_value }))
      if(OPTIMAL_INACCURATE %in% statuses)
        object@status <- OPTIMAL_INACCURATE
      else
        object@status <- OPTIMAL
    }
  }
  object
}

valuesById <- function(object, results_dict, sym_data, solver) {
    if(results_dict[[STATUS]] %in% SOLUTION_PRESENT) {
        outList <- list(value = results_dict[[VALUE]])
        ## Save values for variables in object
        tmp <- saveValuesById(variables(object), sym_data@.var_offsets, results_dict[[PRIMAL]])
        outList <- c(outList, tmp)
        ## Not all solvers provide dual variables
        if(EQ_DUAL %in% names(results_dict)) {
            tmp <- saveDualValues(object, results_dict[[EQ_DUAL]], sym_data@.constr_map[[EQ_MAP]], c("EqConstraint"))
            outList <- c(outList, tmp)
        }
        if(INEQ_DUAL %in% names(results_dict)) {
            tmp <- saveDualValues(object, results_dict[[INEQ_DUAL]], sym_data@.constr_map[[LEQ_MAP]], c("LeqConstraint", "PSDConstraint"))
            outList <- c(outList, tmp)
        }
    } else if(results_dict[[STATUS]] %in% INF_OR_UNB)   # Infeasible or unbounded
        outList <- handleNoSolution(object, results_dict[[STATUS]])
    else   # Solver failed to solve
        stop("Solver failed. Try another.")

    solve_time <- NA_real_
    setup_time <- NA_real_
    num_iters <- NA_real_

    if(SOLVE_TIME %in% names(results_dict))
        solve_time <- results_dict[[SOLVE_TIME]]
    if(SETUP_TIME %in% names(results_dict))
        setup_time <- results_dict[[SETUP_TIME]]
    if(NUM_ITERS %in% names(results_dict))
        num_iters <- results_dict[[NUM_ITERS]]

    result <- c(list(status = results_dict[[STATUS]]),
                outList,
                list(solver = name(solver),
                     solve_time = solve_time,
                     setup_time = setup_time,
                     num_iters = num_iters))
    ##value <- function(cvxObj) result[[ as.character(id(cvxObj)) ]]
    getValue <- function(objet) {
        ## We go French!
        if (is(objet, "Variable") || is(objet, "Constraint")) {
            return(result[[ as.character(id(objet)) ]])
        }
        if(is_zero(objet)) {
            size <- size(objet)
            valResult <- matrix(0, nrow = size[1], ncol = size[2])
        } else {

            arg_values <- list()
            idx <- 1
            for(arg in objet@args) {
                ## An argument without a value makes all higher level values NA.
                ## But if the atom is constant with non-constant arguments, it doesn't depend on its arguments, so it isn't NA.
                arg_val <- if(is_constant(arg)) {
                               value(arg)
                           } else {
                               ##result[[ as.character(id(arg)) ]]
                               getValue(arg)
                           }
                if(is.null(arg_val) || (any(is.na(arg_val)) && !is_constant(objet)))
                    return(NA)
                else {
                    arg_values[[idx]] <- arg_val
                    idx <- idx + 1
                }
            }
            valResult <- to_numeric(objet, arg_values)
        }

        ## Reduce to scalar if possible
        if(all(intf_size(valResult) == c(1, 1)))
            intf_scalar_value(valResult)
        else
            valResult
    }
    getDualValue <- function(objet) {
        if (!is(objet, "Constraint")) {
            stop("getDualValue: argument should be a Constraint!")
        }
        getValue(objet)
    }
    result$getValue <- getValue
    result$getDualValue <- getDualValue
    result
}

# @describeIn Problem Parses the output from a solver and updates the problem state, including the status, objective value, and value of the primal and dual variables. Assumes the results are from the given solver.
# @param object A \linkS4class{Problem} object.
# @param solver A string indicating the name of the solver being used.
# @param results_dict A list containing the solver output.
setMethod("unpack_results", "Problem", function(object, solver, results_dict) {
  canon <- canonicalize(object)
  objective <- canon[[1]]
  constraints <- canon[[2]]
  object@.cached_data <- get_sym_data(solver, objective, constraints, object@.cached_data)
  sym_data <- object@.cached_data[[name(solver)]]@sym_data
  data <- list(sym_data@dims, 0)
  names(data) <- c(DIMS, OFFSET)
  results_dict <- format_results(solver, results_dict, data, object@.cached_data)
  .update_problem_state(object, results_dict, sym_data, solver)
})

handleNoSolution <- function(object, status) {
    ## Set all primal and dual variable values to NA
    ## TODO: This won't work since R creates copies

    result <- list()
    for(var_ in variables(object))
        result[[as.character(id(var_))]] <- NA

    for(constraint in object@constraints)
        result[[as.character(id(constraint))]] <- NA

    ## Set the problem value
    if(tolower(status) %in% c("infeasible", "infeasible_inaccurate"))
        result[[VALUE]] <- primal_to_result(object@objective, Inf)
    else if(tolower(status) %in% c("unbounded", "unbounded_inaccurate"))
        result[[VALUE]] <- primal_to_result(object@objective, -Inf)
    result
}

saveDualValues <- function(object, result_vec, constraints, constr_types) {
    constr_offsets <- integer(0)
    offset <- 0L
    for(constr in constraints) {
        constr_offsets[as.character(constr_id(constr))] <- offset
        offset <- offset + prod(size(constr))
    }

    active_constraints <- list()
    for(constr in object@constraints) {
                                        # Ignore constraints of the wrong type
        if(class(constr) %in% constr_types)
            active_constraints <- c(active_constraints, constr)
    }
    saveValuesById(active_constraints, constr_offsets, result_vec)
}

saveValuesById <- function(variables, offset_map, result_vec) {
    offset_names <- names(offset_map)
    result <- lapply(variables, function(var) {
        id <- as.character(id(var))
        size <- size(var)
        rows <- size[1L]
        cols <- size[2L]
        if (id %in% offset_names) {
            offset <- offset_map[[id]]
            ## Handle scalars
            if (all(c(rows, cols) == c(1L , 1L))) {
                value <- result_vec[offset + 1L]
            } else {
                value <- matrix(result_vec[(offset + 1L):(offset + rows * cols)],
                                nrow = rows, ncol = cols)
            }
        } else {
            ## The variable was multiplied by zero
            value <- matrix(0, nrow = rows, ncol = cols)
        }
        value
    })
    names(result) = sapply(variables, function(x) as.character(id(x)))
    result
 }

## getValue <- function(cvxObj, solution) {
##     solution[[ as.character(id(cvxObj)) ]]
## }

# @describeIn Problem Saves the values of the optimal primal and dual variables.
# @param object A \linkS4class{Problem} object.
# @param result_vec A vector containing the variable values.
# @param objstore The variables or constraints where the values will be stored.
# @param offset_map A list mapping object ID to offset in \code{result_vec}.
setMethod("Problem.save_values", "Problem", function(object, result_vec, objstore, offset_map) {
##  if(length(result_vec) > 0)   # Cast to desired matrix type
##    result_vec <- as.matrix(result_vec)

  objects_new <- list()
  for(obj in objstore) {
    size <- size(obj)
    rows <- size[1]
    cols <- size[2]

    id_char <- as.character(id(obj))
    if(id_char %in% names(offset_map)) {
      offset <- offset_map[[id_char]]

      # Handle scalars
      if(all(c(rows, cols) == c(1,1)))
        # value <- intf_index(result_vec, c(offset, 0))
        value <- result_vec[offset + 1]
      else {
        # value <- matrix(0, nrow = rows, ncol = cols)
        # value <- intf_block_add(value, result_vec[(offset + 1):(offset + rows*cols)], 0, 0, rows, cols)
        value <- result_vec[(offset + 1):(offset + rows*cols)]
        value <- matrix(value, nrow = rows, ncol = cols)
        offset <- offset + rows * cols
      }
    } else  # The variable was multiplied by zero
      value <- matrix(0, nrow = rows, ncol = cols)
    obj <- save_value(obj, value)
    objects_new[[id_char]] <- obj
  }
  objects_new
  # TODO: Update variable values in original problem object
})

#'
#' Arithmetic Operations on Problems
#' 
#' Add, subtract, multiply, or divide DCP optimization problems.
#' 
#' @param e1 The left-hand \linkS4class{Problem} object.
#' @param e2 The right-hand \linkS4class{Problem} object.
#' @return A \linkS4class{Problem} object.
#' @name Problem-arith
NULL

#' @rdname Problem-arith
setMethod("+", signature(e1 = "Problem", e2 = "missing"), function(e1, e2) { Problem(objective = e1@objective, constraints = e1@constraints) })

#' @rdname Problem-arith
setMethod("-", signature(e1 = "Problem", e2 = "missing"), function(e1, e2) { Problem(objective = -e1@objective, constraints = e1@constraints) })

#' @rdname Problem-arith
setMethod("+", signature(e1 = "Problem", e2 = "numeric"), function(e1, e2) { if(length(e2) == 1 && e2 == 0) e1 else stop("Unimplemented") })

#' @rdname Problem-arith
setMethod("+", signature(e1 = "numeric", e2 = "Problem"), function(e1, e2) { e2 + e1 })

#' @rdname Problem-arith
setMethod("+", signature(e1 = "Problem", e2 = "Problem"), function(e1, e2) {
  Problem(objective = e1@objective + e2@objective, constraints = unique(c(e1@constraints, e2@constraints)))
})

#' @rdname Problem-arith
setMethod("-", signature(e1 = "Problem", e2 = "numeric"), function(e1, e2) { e1 + (-e2) })

#' @rdname Problem-arith
setMethod("-", signature(e1 = "numeric", e2 = "Problem"), function(e1, e2) { if(length(e1) == 1 && e2 == 0) -e2 else stop("Unimplemented") })

#' @rdname Problem-arith
setMethod("-", signature(e1 = "Problem", e2 = "Problem"), function(e1, e2) {
  Problem(objective = e1@objective - e2@objective, constraints = unique(c(e1@constraints, e2@constraints)))
})

#' @rdname Problem-arith
setMethod("*", signature(e1 = "Problem", e2 = "numeric"), function(e1, e2) {
  Problem(objective = e1@objective * e2, constraints = e1@constraints)
})

#' @rdname Problem-arith
setMethod("*", signature(e1 = "numeric", e2 = "Problem"), function(e1, e2) { e2 * e1 })

#' @rdname Problem-arith
setMethod("/", signature(e1 = "Problem", e2 = "numeric"), function(e1, e2) {
  Problem(objective = e1@objective * (1.0/e2), constraints = e1@constraints)
})
