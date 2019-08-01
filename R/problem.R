#'
#' The Objective class.
#'
#' This class represents an optimization objective.
#'
#' @slot expr A scalar \linkS4class{Expression} to optimize.
#' @name Objective-class
#' @aliases Objective
#' @rdname Objective-class
.Objective <- setClass("Objective", representation(expr = "ConstValORExpr"), contains = "Canonical")

#' @param expr A scalar \linkS4class{Expression} to optimize.
#' @rdname Objective-class
Objective <- function(expr) { .Objective(expr = expr) }

setMethod("initialize", "Objective", function(.Object, ..., expr) {
  .Object@expr <- expr
  .Object@args <- list(as.Constant(expr))
  # .Object <- callNextMethod(.Object, ..., args = list(as.Constant(expr)))
  
  # Validate that the objective resolves to a scalar.
  if(!is_scalar(.Object@args[[1]]))
    stop("The objective must resolve to a scalar")
  if(!is_real(.Object@args[[1]]))
    stop("The objective must be real valued")
  return(.Object)
})

#' @describeIn Objective The value of the objective expression.
setMethod("value", "Objective", function(object) { 
  v <- value(object@args[[1]])
  if(is.na(v))
    return(NA_real_)
  else
    return(intf_scalar_value(v))
})

#' @describeIn Objective Is the objective a quadratic function?
setMethod("is_quadratic", "Objective", function(object) {
  is_quadratic(object@args[[1]])
})

#' @describeIn Objective Is the objective a quadratic of piecewise affine function?
setMethod("is_qpwa", "Objective", function(object) {
  is_qpwa(object@args[[1]])
})

#'
#' The Minimize class.
#'
#' This class represents an optimization objective for minimization.
#'
#' @slot expr A scalar \linkS4class{Expression} to minimize.
#' @name Minimize-class
#' @aliases Minimize
#' @rdname Minimize-class
.Minimize <- setClass("Minimize", representation(expr = "ConstValORExpr"), contains = "Objective")

#' @param expr A scalar \linkS4class{Expression} to minimize.
#' @rdname Minimize-class
#' @export
Minimize <- function(expr) { .Minimize(expr = expr) }

#' @param object A \linkS4class{Minimize} object.
#' @describeIn Minimize Pass on the target expression's objective and constraints.
setMethod("canonicalize", "Minimize", function(object) { canonical_form(object@args[[1]]) })

#' @describeIn Minimize A logical value indicating whether the objective is convex.
setMethod("is_dcp", "Minimize", function(object) { is_convex(object@args[[1]]) })

#' @describeIn Minimize A logical value indicating whether the objective is log-log convex.
setMethod("is_dgp", "Minimize", function(object) { is_log_log_convex(object@args[[1]]) })

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
#' @name Maximize-class
#' @aliases Maximize
#' @rdname Maximize-class
.Maximize <- setClass("Maximize", contains = "Objective")

#' @param expr A scalar \linkS4class{Expression} to maximize.
#' @rdname Maximize-class
#' @export
Maximize <- function(expr) { .Maximize(expr = expr) }

#' @param object A \linkS4class{Maximize} object.
#' @describeIn Maximize Negates the target expression's objective.
setMethod("canonicalize", "Maximize", function(object) {
  canon <- canonical_form(object@args[[1]])
  list(lo.neg_expr(canon[[1]]), canon[[2]])
})

#' @describeIn Maximize A logical value indicating whether the objective is concave.
setMethod("is_dcp", "Maximize", function(object) { is_concave(object@args[[1]]) })

#' @describeIn Maximize A logical value indicating whether the objective is log-log concave.
setMethod("is_dgp", "Maximize", function(object) { is_log_log_concave(object@args[[1]]) })

# The value of the objective given the solver primal value.
setMethod("primal_to_result", "Maximize", function(object, result) { -result })

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
setMethod("+", signature(e1 = "Objective", e2 = "numeric"), function(e1, e2) { if(length(e2) == 1 && e2 == 0) e1 else stop("Unimplemented") })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "numeric", e2 = "Objective"), function(e1, e2) { e2 + e1 })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Minimize", e2 = "missing"), function(e1, e2) { Maximize(expr = -e1@args[[1]]) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Minimize", e2 = "Minimize"), function(e1, e2) { Minimize(e1@args[[1]] + e2@args[[1]]) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Minimize", e2 = "Maximize"), function(e1, e2) { stop("Problem does not follow DCP rules") })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Objective", e2 = "Minimize"), function(e1, e2) { e1 + (-e2) })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Objective", e2 = "Maximize"), function(e1, e2) { e1 + (-e2) })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Minimize", e2 = "Objective"), function(e1, e2) { (-e2) + e1 })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Maximize", e2 = "Objective"), function(e1, e2) { (-e2) + e1 })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Objective", e2 = "numeric"), function(e1, e2) { if(length(e2) == 1 && e2 == 0) e1 else stop("Unimplemented") })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "numeric", e2 = "Objective"), function(e1, e2) { e2 - e1 })

#' @rdname Objective-arith
setMethod("*", signature(e1 = "Minimize", e2 = "numeric"), function(e1, e2) {
  if(e2 >= 0) Minimize(expr = e1@args[[1]] * e2) else Maximize(expr = e1@args[[1]] * e2)
})

#' @rdname Objective-arith
setMethod("*", signature(e1 = "Maximize", e2 = "numeric"), function(e1, e2) {
  if(e2 < 0) Minimize(expr = e1@args[[1]] * e2) else Maximize(expr = e1@args[[1]] * e2)
})

#' @rdname Objective-arith
setMethod("*", signature(e1 = "numeric", e2 = "Minimize"), function(e1, e2) { e2 * e1 })

#' @rdname Objective-arith
setMethod("*", signature(e1 = "numeric", e2 = "Maximize"), function(e1, e2) { e2 * e1 })

#' @rdname Objective-arith
setMethod("/", signature(e1 = "Objective", e2 = "numeric"), function(e1, e2) { e1 * (1.0/e2) })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Maximize", e2 = "missing"), function(e1, e2) { Minimize(expr = -e1@args[[1]]) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Maximize", e2 = "Maximize"), function(e1, e2) { Maximize(expr = e1@args[[1]] + e2@args[[1]]) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Maximize", e2 = "Minimize"), function(e1, e2) { stop("Problem does not follow DCP rules") })

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
#' @aliases SolverStats
#' @rdname SolverStats-class
.SolverStats <- setClass("SolverStats", representation(solver_name = "character", solve_time = "numeric", setup_time = "numeric", num_iters = "numeric"),
                         prototype(solver_name = NA_character_, solve_time = NA_real_, setup_time = NA_real_, num_iters = NA_real_))

#' @param results_dict A list containing the results returned by the solver.
#' @param solver_name The name of the solver.
#' @return A list containing
#' \describe{
#'   \item{\code{solver_name}}{The name of the solver.}
#'   \item{\code{solve_time}}{The time (in seconds) it took for the solver to solve the problem.}
#'   \item{\code{setup_time}}{The time (in seconds) it took for the solver to set up the problem.}
#'   \item{\code{num_iters}}{The number of iterations the solver had to go through to find a solution.}
#' }
#' @rdname SolverStats-class
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
    
    solver_stats <- list(solver_name, solve_time, setup_time, num_iters)
    names(solver_stats) <- c(SOLVER_NAME, SOLVE_TIME, SETUP_TIME, NUM_ITERS)
    return(solver_stats)
    ## .SolverStats(solver_name = solver_name, solve_time = solve_time, setup_time = setup_time, num_iters = num_iters)
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
#' @aliases SizeMetrics
#' @rdname SizeMetrics-class
.SizeMetrics <- setClass("SizeMetrics", representation(num_scalar_variables = "numeric", num_scalar_data = "numeric", num_scalar_eq_constr = "numeric", num_scalar_leq_constr = "numeric",
                                                       max_data_dimension = "numeric", max_big_small_squared = "numeric"),
                         prototype(num_scalar_variables = NA_real_, num_scalar_data = NA_real_, num_scalar_eq_constr = NA_real_, num_scalar_leq_constr = NA_real_,
                                   max_data_dimension = NA_real_, max_big_small_squared = NA_real_))

#' @param problem A \linkS4class{Problem} object.
#' @rdname SizeMetrics-class
SizeMetrics <- function(problem) {
  # num_scalar_variables
  num_scalar_variables <- 0
  for(var in variables(problem))
    num_scalar_variables <- num_scalar_variables + size(var)

  # num_scalar_data, max_data_dimension, and max_big_small_squared
  max_data_dimension <- 0
  num_scalar_data <- 0
  max_big_small_squared <- 0
  for(const in c(constants(problem), parameters(problem))) {
    big <- 0
    # Compute number of data
    num_scalar_data <- num_scalar_data + size(const)
    big <- ifelse(1, length(dim(const)) == 0, max(dim(const)))
    small <- ifelse(1, length(dim(const)) == 0, min(dim(const)))

    # Get max data dimension
    if(max_data_dimension < big)
      max_data_dimension <- big

    if(max_big_small_squared < as.numeric(big)*small*small)
      max_big_small_squared <- as.numeric(big)*small*small
  }

  # num_scalar_eq_constr
  num_scalar_eq_constr <- 0
  for(constraint in problem@constraints) {
    if(is(constraint, "EqConstraint") || is(constraint, "ZeroConstraint"))
      num_scalar_eq_constr <- num_scalar_eq_constr + size(constraint@expr)
  }

  # num_scalar_leq_constr
  num_scalar_leq_constr <- 0
  for(constraint in problem@constraints) {
    if(is(constraint, "IneqConstraint") || is(constraint, "NonPosConstraint"))
      num_scalar_leq_constr <- num_scalar_leq_constr + size(constraint@expr)
  }

  .SizeMetrics(num_scalar_variables = num_scalar_variables, num_scalar_data = num_scalar_data, num_scalar_eq_constr = num_scalar_eq_constr, num_scalar_leq_constr = num_scalar_leq_constr,
               max_data_dimension = max_data_dimension, max_big_small_squared = max_big_small_squared)
}

setClassUnion("SizeMetricsORNULL", c("SizeMetrics", "NULL"))

#'
#' The Solution class.
#'
#' This class represents a solution to an optimization problem.
#'
#' @rdname Solution-class
.Solution <- setClass("Solution", representation(status = "character", opt_val = "numeric", primal_vars = "list", dual_vars = "list", attr = "list"),
                      prototype(primal_vars = list(), dual_vars = list(), attr = list()))

Solution <- function(status, opt_val, primal_vars, dual_vars, attr) {
  .Solution(status = status, opt_val = opt_val, primal_vars = primal_vars, dual_vars = dual_vars, attr = attr)
}

setMethod("show", "Solution", function(object) {
  cat("Solution(", object@status, ", (",
      paste(object@primal_vars, collapse = ", "), "), (",
      paste(object@dual_vars, collapse = ", "), "), (",
      paste(object@attr, collapse = ", "), "))", sep = "")
})

# TODO: Get rid of this and just skip calling copy on Solution objects.
setMethod("copy", "Solution", function(object, args = NULL, id_objects = list()) { return(object) })

setMethod("as.character", "Solution", function(x) {
  paste("Solution(", x@status, ", (",
        paste(x@primal_vars, collapse = ", "), "), (",
        paste(x@dual_vars, collapse = ", "), "), (",
        paste(x@attr, collapse = ", "), "))", sep = "")
})

setClassUnion("SolutionORList", c("Solution", "list"))

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
#' @aliases Problem
#' @rdname Problem-class
.Problem <- setClass("Problem", representation(objective = "Objective", constraints = "list", variables = "list", value = "numeric", status = "character", solution = "ANY", .intermediate_chain = "ANY", .solving_chain = "ANY", .cached_chain_key = "list", .separable_problems = "list", .size_metrics = "SizeMetricsORNULL", .solver_stats = "list", args = "list", .solver_cache = "list", .intermediate_problem = "ANY", .intermediate_inverse_data = "ANY"),
                    prototype(constraints = list(), value = NA_real_, status = NA_character_, solution = NULL, .intermediate_chain = NULL, .solving_chain = NULL, .cached_chain_key = list(), .separable_problems = list(), .size_metrics = NULL, .solver_stats = NULL, args = list(), .solver_cache = list(), .intermediate_problem = NULL, .intermediate_inverse_data = NULL),
                    validity = function(object) {
                      if(!(class(object@objective) %in% c("Minimize", "Maximize")))
                        stop("[Problem: objective] objective must be Minimize or Maximize")
                      if(!is.na(object@value))
                        stop("[Problem: value] value should not be set by user")
                      if(!is.na(object@status))
                        stop("[Problem: status] status should not be set by user")
                      if(!is.null(object@solution))
                        stop("[Problem: solution] solution should not be set by user")
                      if(!is.null(object@.intermediate_chain))
                        stop("[Problem: .intermediate_chain] .intermediate_chain is an internal slot and should not be set by user")
                      if(!is.null(object@.solving_chain))
                        stop("[Problem: .solving_chain] .solving_chain is an internal slot and should not be set by user")
                      if(length(object@.cached_chain_key) > 0)
                        stop("[Problem: .cached_chain_key] .cached_chain_key is an internal slot and should not be set by user")
                      if(length(object@.separable_problems) > 0)
                        stop("[Problem: .separable_problems] .separable_problems is an internal slot and should not be set by user")
                      if(!is.null(object@.size_metrics))
                        stop("[Problem: .size_metrics] .size_metrics is an internal slot and should not be set by user")
                      if(length(object@.solver_stats) > 0)
                        stop("[Problem: .solver_stats] .solver_stats is an internal slot and should not be set by user")
                      if(length(object@args) > 0)
                        stop("[Problem: args] args is an internal slot and should not be set by user")
                      if(length(object@.solver_cache) > 0)
                        stop("[Problem: .solver_cache] .solver_cache is an internal slot and should not be set by user")
                      if(!is.null(object@.intermediate_problem))
                        stop("[Problem: .intermediate_problem] .intermediate_problem is an internal slot and should not be set by user")
                      if(!is.null(object@.intermediate_inverse_data))
                        stop("[Problem: .intermediate_inverse_data] .intermediate_inverse_data is an internal slot and should not be set by user")
                      return(TRUE)
                    }, contains = "Canonical")

#' @param objective A \linkS4class{Minimize} or \linkS4class{Maximize} object representing the optimization objective.
#' @param constraints (Optional) A list of \linkS4class{Constraint} objects representing constraints on the optimization variables.
#' @rdname Problem-class
#' @export
Problem <- function(objective, constraints = list()) {
  .Problem(objective = objective, constraints = constraints)
}

# Used by pool.map to send solve result back. Unsure if this is necessary for multithreaded operation in R.
SolveResult <- function(opt_value, status, primal_values, dual_values) { list(opt_value = opt_value, status = status, primal_values = primal_values, dual_values = dual_values, class = "SolveResult") }

setMethod("initialize", "Problem", function(.Object, ..., objective, constraints = list(), variables, value = NA_real_, status = NA_character_, solution = NULL, .intermediate_chain = NULL, .solving_chain = NULL, .cached_chain_key = list(), .separable_problems = list(), .size_metrics = SizeMetrics(), .solver_stats = list(), args = list(), .solver_cache = list(), .intermediate_problem = NULL, .intermediate_inverse_data = NULL) {
  .Object@objective <- objective
  .Object@constraints <- constraints
  .Object@variables <- Problem.build_variables(.Object)
  .Object@value <- value
  .Object@status <- status
  .Object@solution <- solution
  
  .Object@.intermediate_problem <- .intermediate_problem
  .Object@.intermediate_inverse_data <- .intermediate_inverse_data
  
  # The intermediate and solving chains to canonicalize and solve the problem.
  .Object@.intermediate_chain <- .intermediate_chain
  .Object@.solving_chain <- .solving_chain
  .Object@.cached_chain_key <- .cached_chain_key

  # List of separable (sub)problems
  .Object@.separable_problems <- .separable_problems

  # Information about the dimensions of the problem and its constituent parts
  .Object@.size_metrics <- SizeMetrics(.Object)

  # Benchmarks reported by the solver.
  .Object@.solver_stats <- .solver_stats
  .Object@args <- list(.Object@objective, .Object@constraints)
  
  # Cache for warm start.
  .Object@.solver_cache <- list()
  .Object
})

#' @param object A \linkS4class{Problem} object.
#' @describeIn Problem The objective of the problem.
setMethod("objective", "Problem", function(object) { object@objective })

#' @describeIn Problem Set the value of the problem objective.
setReplaceMethod("objective", "Problem", function(object, value) {
  object@objective <- value
  object
})

#' @describeIn Problem A list of the constraints of the problem.
setMethod("constraints", "Problem", function(object) { object@constraints })

#' @describeIn Problem Set the value of the problem constraints.
setReplaceMethod("constraints", "Problem", function(object, value) {
  object@constraints <- value
  object
})

#' @describeIn Problem The value from the last time the problem was solved (or NA if not solved).
setMethod("value", "Problem", function(object) { 
  if(is.na(object@value))
    return(NA_real_)
  else
    return(intf_scalar_value(object@value))
})

#' @param value A \linkS4class{Minimize} or \linkS4class{Maximize} object (objective), list of \linkS4class{Constraint} objects (constraints), or numeric scalar (value).
#' @describeIn Problem Set the value of the optimal objective.
setReplaceMethod("value", "Problem", function(object, value) {
    object@value <- value
    object
})

#' @describeIn Problem The status from the last time the problem was solved.
setMethod("status", "Problem", function(object) { object@status })

# Set the status of the problem.
setReplaceMethod("status", "Problem", function(object, value) {
    object@status <- value
    object
})

#' @describeIn Problem A logical value indicating whether the problem statisfies DCP rules.
#' @examples
#' x <- Variable(2)
#' p <- Problem(Minimize(p_norm(x, 2)), list(x >= 0))
#' is_dcp(p)
setMethod("is_dcp", "Problem", function(object) {
  all(sapply(c(object@constraints, list(object@objective)), is_dcp))
})

#' @describeIn Problem A logical value indicating whether the problem statisfies DGP rules.
setMethod("is_dgp", "Problem", function(object) {
  all(sapply(c(object@constraints, list(object@objective)), is_dgp))
})

#' @describeIn Problem A logical value indicating whether the problem is a quadratic program.
#' @examples
#' x <- Variable(2)
#' A <- matrix(c(1,-1,-1, 1), nrow = 2)
#' p <- Problem(Minimize(quad_form(x, A)), list(x >= 0))
#' is_qp(p)
setMethod("is_qp", "Problem", function(object) {
  for(c in object@constraints) {
    if(!(is(c, "EqConstraint") || is(c, "ZeroConstraint") || is_pwl(c@args[[1]])))
      return(FALSE)
  }
  for(var in variables(object)) {
    if(is_psd(var) || is_nsd(var))
      return(FALSE)
  }
  return(is_dcp(object) && is_qpwa(object@objective@args[[1]]))
})

#' @describeIn Problem The graph implementation of the problem.
setMethod("canonicalize", "Problem", function(object) {
  obj_canon <- canonical_form(object@objective)
  canon_constr <- obj_canon[[2]]

  for(constr in object@constraints)
    canon_constr <- c(canon_constr, canonical_form(constr)[[2]])
  list(obj_canon[[1]], canon_constr)
})

setMethod("is_mixed_integer", "Problem", function(object) {
  any(sapply(variables(object), function(v) { v@attributes$boolean || v@attributes$integer }))
})

#' @describeIn Problem List of \linkS4class{Variable} objects in the problem.
setMethod("variables", "Problem", function(object) { object@variables })

Problem.build_variables <- function(object) {
  vars_ <- variables(object@objective)
  constrs_ <- lapply(object@constraints, function(constr) { variables(constr) })
  unique(flatten_list(c(vars_, constrs_)))   # Remove duplicates
}

#' @describeIn Problem List of \linkS4class{Parameter} objects in the problem.
setMethod("parameters", "Problem", function(object) {
  params <- parameters(object@objective)
  constrs_ <- lapply(object@constraints, function(constr) { parameters(constr) })
  unique(flatten_list(c(params, constrs_)))   # Remove duplicates
})

#' @describeIn Problem List of \linkS4class{Constant} objects in the problem.
setMethod("constants", "Problem", function(object) {
  constants_ <- lapply(object@constraints, function(constr) { constants(constr) })
  constants_ <- c(constants(object@objective), constants_)
  unique(flatten_list(constants_))   # TODO: Check duplicated constants are removed correctly
})

#' @describeIn Problem List of \linkS4class{Atom} objects in the problem.
setMethod("atoms", "Problem", function(object) {
  atoms_ <- lapply(object@constraints, function(constr) { atoms(constr) })
  atoms_ <- c(atoms_, atoms(object@objective))
  unique(flatten_list(atoms_))
})

#' @describeIn Problem Information about the size of the problem.
setMethod("size_metrics", "Problem", function(object) { object@.size_metrics })

# Additional information returned by the solver.
setMethod("solver_stats", "Problem", function(object) { object@.solver_stats })

# Set the additional information returned by the solver in the problem.
setMethod("solver_stats<-", "Problem", function(object, value) {
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

#' @param solver A string indicating the solver that the problem data is for. Call \code{installed_solvers()} to see all available.
#' @describeIn Problem Get the problem data passed to the specified solver.
setMethod("get_problem_data", signature(object = "Problem", solver = "character", gp = "logical"), function(object, solver, gp) {
  object <- .construct_chains(object, solver = solver, gp = gp)
  
  tmp <- perform(object@.solving_chain, object@.intermediate_problem)
  data <- tmp[[1]]
  solving_inverse_data <- tmp[[2]]
  
  full_chain <- prepend(object@.solving_chain, object@.intermediate_chain)
  inverse_data <- c(object@.intermediate_inverse_data, solving_inverse_data)
  
  return(list(data, full_chain, inverse_data))
})

.find_candidate_solvers <- function(object, solver = NA, gp = FALSE) {
  candidates <- list(qp_solvers = list(), conic_solvers = list())
  
  if(!is.na(solver)) {
    if(!(solver %in% INSTALLED_SOLVERS))
      stop("The solver is not installed")
    if(solver %in% CONIC_SOLVERS)
      candidates$conic_solvers <- c(candidates$conic_solvers, solver)
    if(solver %in% QP_SOLVERS)
      candidates$qp_solvers <- c(candidates$qp_solvers, solver)
  } else {
    candidates$qp_solvers <- INSTALLED_SOLVERS[INSTALLED_SOLVERS %in% QP_SOLVERS]
    candidates$conic_solvers <- INSTALLED_SOLVERS[INSTALLED_SOLVERS %in% CONIC_SOLVERS]
  }
  
  # If gp, we must have only conic solvers.
  if(gp) {
    if(!is.na(solver) && !(solver %in% CONIC_SOLVERS))
      stop("When gp = TRUE, solver must be a conic solver. Try calling solve() with solver = ECOS()")
    else if(is.na(solver))
      candidates$qp_solvers <- list()   # No QP solvers allowed.
  }
  
  if(is_mixed_integer(object)) {
    if(length(candidates$qp_solvers) > 0) {
      qp_filter <- sapply(candidates$qp_solvers, function(s) { mip_capable(SOLVER_MAP_QP[[s]]) })
      candidates$qp_solvers <- candidates$qp_solvers[qp_filter]
    }
    
    if(length(candidates$conic_solvers) > 0) {
      conic_filter <- sapply(candidates$conic_solvers, function(s) { mip_capable(SOLVER_MAP_CONIC[[s]]) })
      candidates$conic_solvers <- candidates$conic_solvers[conic_filter]
    }
    
    if(length(candidates$conic_solvers) == 0 && length(candidates$qp_solvers) == 0)
      stop("Problem is mixed-integer, but candidate QP/Conic solvers are not MIP-capable")
  }
  return(candidates)
}

setMethod("get_problem_data", signature(object = "Problem", solver = "character", gp = "missing"), function(object, solver, gp) {
  get_problem_data(object, solver, gp = FALSE)
})

.construct_chains <- function(object, solver = NA, gp = FALSE) {
  chain_key <- list(solver, gp)
  
  if(!identical(chain_key, object@.cached_chain_key)) {
    candidate_solvers <- .find_candidate_solvers(object, solver = solver, gp = gp)
    object@.intermediate_chain <- construct_intermediate_chain(object, candidate_solvers, gp = gp)
    tmp <- perform(object@.intermediate_chain, object)
    object@.intermediate_problem <- tmp[[1]]
    object@.intermediate_inverse_data <- tmp[[2]]
    
    object@.solving_chain <- construct_solving_chain(object@.intermediate_problem, candidate_solvers)
    object@.cached_chain_key <- chain_key
  }
  return(object)
}

#' @docType methods
#' @rdname psolve
#' @export
setMethod("psolve", "Problem", function(object, solver = NA, ignore_dcp = FALSE, warm_start = FALSE, verbose = FALSE, parallel = FALSE, gp = FALSE, ...) {
  if(parallel)
    stop("Unimplemented")
  
  object <- .construct_chains(object, solver = solver, gp = gp)
  tmp <- perform(object@.solving_chain, object@.intermediate_problem)
  data <- tmp[[1]]
  solving_inverse_data <- tmp[[2]]
  solution <- reduction_solve_via_data(object@.solving_chain, object, data, warm_start, verbose, list(...))
  
  full_chain <- prepend(object@.solving_chain, object@.intermediate_chain)
  inverse_data <- c(object@.intermediate_inverse_data, solving_inverse_data)
  # object <- unpack_results(object, solution, full_chain, inverse_data)
  # return(value(object))
  unpack_results(object, solution, full_chain, inverse_data)
})

#' @docType methods
#' @rdname psolve
#' @method solve Problem
#' @export
setMethod("solve", signature(a = "Problem", b = "ANY"), function(a, b = NA, ...) {
  kwargs <- list(...)
  if(missing(b)) {
    if("solver" %in% names(kwargs))
      psolve(a, ...)
    else
      psolve(a, solver = NA, ...)
  } else
    psolve(a, b, ...)
})

# # TODO: Finish implementation of parallel solve.
# .parallel_solve.Problem <- function(object, solver = NULL, ignore_dcp = FALSE, warm_start = FALSE, verbose = FALSE, ...) {
#   .solve_problem <- function(problem) {
#     result <- solve(problem, solver = solver, ignore_dcp = ignore_dcp, warm_start = warm_start, verbose = verbose, parallel = FALSE, ...)
#     opt_value <- result$optimal_value
#     status <- result$status
#     primal_values <- result$primal_values
#     dual_values <- result$dual_values
#     return(SolveResults(opt_value, status, primal_values, dual_values))
#   }
#
#   # TODO: Finish implementation of parallel solve
#   statuses <- sapply(solve_results, function(solve_result) { solve_result$status })
#
#   # Check if at least one subproblem is infeasible or inaccurate
#   for(status in INF_OR_UNB) {
#     if(status %in% statuses) {
#       .handle_no_solution(object, status)
#       break
#     } else {
#       for(i in 1:length(solve_results)) {
#         subproblem <- object@separable_problems[i]
#         solve_result <- solve_results[i]
#
#         for(j in 1:length(solve_result$primal_values)) {
#           var <- variables(subproblem)[j]
#           primal_value <- solve_result$primal_values[j]
#           subproblem <- save_value(var, primal_value)   # TODO: Fix this since R makes copies
#         }
#
#         for(j in 1:length(solve_result$dual_values)) {
#           constr <- subproblem@constraints[j]
#           dual_value <- solve_result$dual_values[j]
#           subproblem <- save_value(constr, dual_value)   # TODO: Fix this since R makes copies
#         }
#       }
#
#       object@value <- sum(sapply(solve_results, function(solve_result) { solve_result$optimal_value }))
#       if(OPTIMAL_INACCURATE %in% statuses)
#         object@status <- OPTIMAL_INACCURATE
#       else
#         object@status <- OPTIMAL
#     }
#   }
#   object
# }

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
      if(is(objet, "Variable") || is(objet, "Constraint"))
        return(result[[as.character(id(objet))]])
      if(is_zero(objet)) {
        dims <- dim(objet)
        valResult <- matrix(0, nrow = dims[1], ncol = dims[2])
      } else {
        arg_values <- list()
        idx <- 1
        for(arg in objet@args) {
          ## An argument without a value makes all higher level values NA.
          ## But if the atom is constant with non-constant arguments, it doesn't depend on its arguments, so it isn't NA.
          arg_val <- if(is_constant(arg))
                       value(arg)
                     else {
                       ## result[[as.character(id(arg))]]
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
      if(all(intf_dim(valResult) == c(1, 1)))
        intf_scalar_value(valResult)
      else
        valResult
    }
    
    getDualValue <- function(objet) {
      if(!is(objet, "Constraint")) {
        stop("getDualValue: argument should be a Constraint!")
      }
      getValue(objet)
    }
    
    result$getValue <- getValue
    result$getDualValue <- getDualValue
    result
}

# .clear_solution.Problem <- function(object) {
#   for(v in variables(object))
#     object <- save_value(v, NA)
#   for(c in constraints(object))
#     object <- save_value(c, NA)
#   object@value <- NA
#   object@status <- NA
#   object@solution <- NA
#   return(object)
# }

setMethod("unpack", signature(object = "Problem", solution = "Solution"), function(object, solution) {
  # if(solution@status %in% SOLUTION_PRESENT) {
  #   for(v in variables(object))
  #     object <- save_value(v, solution@primal_vars[id(v)])
  #   for(c in object@constraints) {
  #     if(id(c) %in% solution@dual_vars)
  #       object <- save_value(c, solution@dual_vars[id(c)])
  #   }
  # } else if(solution@status %in% INF_OR_UNB) {
  #   for(v in variables(object))
  #     object <- save_value(v, NA)
  #   for(constr in object@constraints)
  #     object <- save_value(constr, NA)
  # } else
  #   stop("Cannot unpack invalid solution")
  # object@value <- solution@opt_val
  # object@status <- solution@status
  # object@solution <- solution
  # return(object)
  
  result <- list()
  if(solution@status %in% SOLUTION_PRESENT) {
    for(v in variables(object)) {
      vid <- as.character(id(v))
      val <- solution@primal_vars[[vid]]
      if(is.null(dim(val)) || all(dim(val) == 1))
        val <- as.vector(val)
      result[[vid]] <- val
      # result[[vid]] <- solution@primal_vars[[vid]]
    }
    for(c in object@constraints) {
      cid <- as.character(id(c))
      if(cid %in% names(solution@dual_vars)) {
        val <- solution@dual_vars[[cid]]
        if(is.null(dim(val)) || all(dim(val)) == 1)
          val <- as.vector(val)
        result[[cid]] <- val
        # result[[cid]] <- solution@dual_vars[[cid]]
      }
    }
  } else if(solution@status %in% INF_OR_UNB) {
    for(v in variables(object))
      result[[as.character(id(v))]] <- NA
    for(c in object@constraints)
      result[[as.character(id(c))]] <- NA
  } else if(solution@status %in% ERROR) {
    warning("Solver returned with status ", solution@status)
    result$status <- solution@status
    return(result)
  } else
    stop("Cannot unpack invalid solution")
  
  result$value <- solution@opt_val
  result$status <- solution@status
  # result$solution <- solution
  
  # Helper functions.
  getValue <- function(objet) {
    ## We go French!
    if(is(objet, "Variable") || is(objet, "Constraint"))
      return(result[[as.character(id(objet))]])
    if(is_zero(objet)) {
      dims <- dim(objet)
      valResult <- matrix(0, nrow = dims[1], ncol = dims[2])
    } else {
      arg_values <- list()
      idx <- 1
      for(arg in objet@args) {
        ## An argument without a value makes all higher level values NA.
        ## But if the atom is constant with non-constant arguments, it doesn't depend on its arguments, so it isn't NA.
        arg_val <- if(is_constant(arg))
          value(arg)
        else {
          ## result[[as.character(id(arg))]]
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
    if(all(intf_dim(valResult) == c(1, 1)))
      intf_scalar_value(valResult)
    else
      valResult
  }
  
  getDualValue <- function(objet) {
    if(!is(objet, "Constraint"))
      stop("getDualValue: argument should be a Constraint!")
    getValue(objet)
  }
  
  result$getValue <- getValue
  result$getDualValue <- getDualValue
  return(result)
})

#' @describeIn Problem
#' Parses the output from a solver and updates the problem state, including the status,
#' objective value, and values of the primal and dual variables.
#' Assumes the results are from the given solver.
#' @param results_dict A list containing the solver output.
#' @docType methods
setMethod("unpack_results", "Problem", function(object, solution, chain, inverse_data) {
  solution <- invert(chain, solution, inverse_data)
  # object <- unpack(object, solution)
  # object@.solver_stats <- SolverStats(object@solution@attr, name(chain@solver))
  # return(object)
  results <- unpack(object, solution)
  solver_stats <- SolverStats(solution@attr, name(chain@solver))
  return(c(results, solver_stats))
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
        constr_offsets[as.character(id(constr))] <- offset
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
