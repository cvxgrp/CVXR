setClassUnion("ListORNULL", c("list", "NULL"))


#'
#' The Objective class.
#'
#' This class represents an optimization objective.
#'
#' @slot expr A scalar \linkS4class{Expression} to optimize.
#' @slot NAME (Internal) Name used for debug messages.
#' @name Objective-class
#' @aliases Objective
#' @rdname Objective-class
.Objective <- setClass("Objective",
                       representation(expr = "ConstValORExpr", NAME = "character"),
                       contains = "Canonical")

#' @param expr A scalar \linkS4class{Expression} to optimize.
#' @rdname Objective-class
Objective <- function(expr) { .Objective(expr = expr) }

setMethod("initialize", "Objective", function(.Object, ..., expr) {
  .Object@expr <- expr
  .Object@NAME  <- "objective"
  .Object@args <- list(as.Constant(expr))
  # .Object <- callNextMethod(.Object, ..., args = list(as.Constant(expr)))

  # Validate that the objective resolves to a scalar.
  if(!is_scalar(.Object@args[[1]]))
    stop("The objective must resolve to a scalar")
  if(!is_real(.Object@args[[1]]))
    stop("The objective must be real valued")
  return(.Object)
})

setMethod("show", "Objective", function(object) {
  cat("Objective(", name(object@args[[1]]), ")", sep = "")
})

setMethod("as.character", "Objective", function(x) {
  paste("Objective", name(x@args[[1]]))
})

#' @param object An \linkS4class{Objective} object.
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
setMethod("is_dcp", "Minimize", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(is_convex(object@args[[1]]))
  }
  return(is_convex(object@args[[1]]))
})

#' @describeIn Minimize A logical value indicating whether the objective is log-log convex.
setMethod("is_dgp", "Minimize", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(is_log_log_convex(object@args[[1]]))
  }
  return(is_log_log_convex(object@args[[1]]))
})

#' @describeIn Minimize Is the objective DPP?
setMethod("is_dpp", "Minimize", function(object, context = "dcp") {
  dpp_scope()   # TODO: Implement DPP scoping.
  if(tolower(context) == "dcp")
    return(is_dcp(object, dpp = TRUE))
  else if(tolower(context) == "dgp")
    return(is_dgp(object, dpp = TRUE))
  else
    stop("Unsupported context ", context)
})

#' @describeIn Minimize A logical value indicating whether the objective is quasiconvex.
setMethod("is_dqcp", "Minimize", function(object) {
  is_quasiconvex(object@args[[1]])
})

# The value of the objective given the solver primal value.
Minimize.primal_to_result <- function(result) { result }

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
setMethod("is_dcp", "Maximize", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(is_concave(object@args[[1]]))
  }
  return(is_concave(object@args[[1]]))
})

#' @describeIn Maximize A logical value indicating whether the objective is log-log concave.
setMethod("is_dgp", "Maximize", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(is_log_log_concave(object@args[[1]]))
  }
  return(is_log_log_concave(object@args[[1]]))
})

#' @describeIn Maximize Is the objective DPP?
setMethod("is_dpp", "Maximize", function(object, context = "dcp") {
  dpp_scope()   # TODO: Implement DPP scoping.
  if(tolower(context) == "dcp")
    return(is_dcp(object, dpp = TRUE))
  else if(tolower(context) == "dgp")
    return(is_dgp(object, dpp = TRUE))
  else
    stop("Unsupported context ", context)
})

#' @describeIn Maximize A logical value indicating whether the objective is quasiconcave.
setMethod("is_dqcp", "Maximize", function(object) { is_quasiconcave(object@args[[1]]) })

# The value of the objective given the solver primal value.
Maximize.primal_to_result <- function(result) { -result }

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

#============================#
# Objective Class Arithmetic
#============================#
#' @rdname Objective-arith
setMethod("+", signature(e1 = "Objective", e2 = "numeric"), function(e1, e2) { if(length(e2) == 1 && e2 == 0) e1 else stop("Unimplemented") })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "numeric", e2 = "Objective"), function(e1, e2) { e2 + e1 })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Objective", e2 = "numeric"), function(e1, e2) { if(length(e2) == 1 && e2 == 0) e1 else stop("Unimplemented") })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "numeric", e2 = "Objective"), function(e1, e2) { if(length(e1) == 0 && e1 == 0) -e2 else stop("Unimplemented") })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Objective", e2 = "Minimize"), function(e1, e2) { e1 + (-e2) })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Objective", e2 = "Maximize"), function(e1, e2) { e1 + (-e2) })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Minimize", e2 = "Objective"), function(e1, e2) { (-e2) + e1 })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Maximize", e2 = "Objective"), function(e1, e2) { (-e2) + e1 })

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

#===========================#
# Minimize Class Arithmetic
#===========================#
#' @rdname Objective-arith
setMethod("-", signature(e1 = "Minimize", e2 = "missing"), function(e1, e2) { Maximize(expr = -e1@args[[1]]) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Minimize", e2 = "Minimize"), function(e1, e2) { Minimize(e1@args[[1]] + e2@args[[1]]) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Minimize", e2 = "Maximize"), function(e1, e2) { stop("Problem does not follow DCP rules") })

#===========================#
# Maximize Class Arithmetic
#===========================#
#' @rdname Objective-arith
setMethod("-", signature(e1 = "Maximize", e2 = "missing"), function(e1, e2) { Minimize(expr = -e1@args[[1]]) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Maximize", e2 = "Maximize"), function(e1, e2) { Maximize(expr = e1@args[[1]] + e2@args[[1]]) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Maximize", e2 = "Minimize"), function(e1, e2) { stop("Problem does not follow DCP rules") })

#setClassUnion("SolvingChainORNULL", c("SolvingChain", "NULL"))
## We add SolvingChain later to the class union below using setIs function
## see reduction_solvers.R near SolvingChain class defn
setClassUnion("SolvingChainORNULL", c("NULL"))

#setClassUnion("ParamProbORNULL", c("ParamProg", "NULL")) ## THIS SEEMS WRONG, perhaps meant next line
#setClassUnion("ParamProgORNULL", c("ParamQuadProg", "ParamConeProg", "NULL"))
## We add ParamQuadProg, ParamConeProg later to class union below using setIs
## see dcp2cone.R near ParamConeProg class defn and dcp2quad.R near ParamQuadProg
setClassUnion("ParamProgORNULL", c("NULL"))

#setClassUnion("InverseDataORNULL", c("InverseData", "NULL"))
## We add InverseData to InverseDataORNULL later to class union below using setIs
## see reductions.R near InverseData class defn
setClassUnion("InverseDataORNULL", c("NULL"))

#'
#' The Cache class.
#'
#' This class contains cached information.
#'
#' @rdname Cache-class
Cache <- setClass("Cache", representation(key = "ListORNULL", solving_chain = "SolvingChainORNULL", param_prog = "ParamProgORNULL", inverse_data = "InverseDataORNULL"),
                  prototype(key = NULL, solving_chain = NULL, param_prog = NULL, inverse_data = NULL))

setMethod("invalidate", "Cache", function(object) {
  object@key <- list()
  object@solving_chain <- NULL
  object@param_prog <- NULL
  object@inverse_data <- NULL
  return(object)
})

setMethod("make_key", "Cache", function(object, solver, gp, ignore_dpp) {
  return(list(solver, gp, ignore_dpp))
})

setMethod("gp", "Cache", function(object) {
  return(!is.null(object@key) && !is.na(object@key) && object@key[[2]])
})

#'
#' The SolverStats class.
#'
#' This class contains the miscellaneous information that is returned by a solver after solving, but that is not captured directly by the \linkS4class{Problem} object.
#'
#' @slot solver_name The name of the solver.
#' @slot solve_time The time (in seconds) it took for the solver to solve the problem.
#' @slot setup_time The time (in seconds) it took for the solver to set up the problem.
#' @slot num_iters The number of iterations the solver had to go through to find a solution.
#' @slot extra_stats Extra statistics specific to the solver; these statistics are typically returned directly from the solver, without modification by CVXR.
#' @name SolverStats-class
#' @aliases SolverStats
#' @rdname SolverStats-class
.SolverStats <- setClass("SolverStats", representation(solver_name = "character", solve_time = "numeric", setup_time = "numeric", num_iters = "numeric", extra_stats = "ListORNULL"),
                         prototype(solver_name = NA_character_, solve_time = NA_real_, setup_time = NA_real_, num_iters = NA_real_, extra_stats = NULL))

#' @param results_dict A list containing the results returned by the solver.
#' @param solver_name The name of the solver.
#' @return A list containing
#' \describe{
#'   \item{\code{solver_name}}{The name of the solver.}
#'   \item{\code{solve_time}}{The time (in seconds) it took for the solver to solve the problem.}
#'   \item{\code{setup_time}}{The time (in seconds) it took for the solver to set up the problem.}
#'   \item{\code{num_iters}}{The number of iterations the solver had to go through to find a solution.}
#'   \item{\code{extra_stats}}{Extra statistics specific to the solver; these statistics are typically returned directly from the solver, without modification by CVXR.}
#' }
#' @rdname SolverStats-class
SolverStats <- function(results_dict = list(), solver_name = NA_character_) {
    solve_time <- NA_real_
    setup_time <- NA_real_
    num_iters <- NA_real_
    extra_stats <- NULL

    if(SOLVE_TIME %in% names(results_dict))
        solve_time <- results_dict[[SOLVE_TIME]]
    if(SETUP_TIME %in% names(results_dict))
        setup_time <- results_dict[[SETUP_TIME]]
    if(NUM_ITERS %in% names(results_dict))
        num_iters <- results_dict[[NUM_ITERS]]
    if(EXTRA_STATS %in% names(results_dict))
      extra_stats <- results_dict[[EXTRA_STATS]]

    .SolverStats(solver_name = solver_name, solve_time = solve_time, setup_time = setup_time, num_iters = num_iters, extra_stats = extra_stats)
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
    big <- ifelse(length(dim(const)) == 0, 1, max(dim(const)))
    small <- ifelse(length(dim(const)) == 0, 1, min(dim(const)))

    # Get max data dimension
    if(max_data_dimension < big)
      max_data_dimension <- big

    max_big_small_squared_init <- as.numeric(big)*as.numeric(small)^2
    if(max_big_small_squared < max_big_small_squared_init)
      max_big_small_squared <- max_big_small_squared_init
  }

  # num_scalar_eq_constr
  num_scalar_eq_constr <- 0
  for(constraint in problem@constraints) {
    if(is(constraint, "EqConstraint") || is(constraint, "ZeroConstraint"))
      num_scalar_eq_constr <- num_scalar_eq_constr + size(expr(constraint))
  }

  # num_scalar_leq_constr
  num_scalar_leq_constr <- 0
  for(constraint in problem@constraints) {
    if(is(constraint, "IneqConstraint") || is(constraint, "NonPosConstraint") || is(constraint, "NonNegConstraint"))
      num_scalar_leq_constr <- num_scalar_leq_constr + size(expr(constraint))
  }

  .SizeMetrics(num_scalar_variables = num_scalar_variables, num_scalar_data = num_scalar_data, num_scalar_eq_constr = num_scalar_eq_constr, num_scalar_leq_constr = num_scalar_leq_constr,
               max_data_dimension = max_data_dimension, max_big_small_squared = max_big_small_squared)
}

setClassUnion("SolverStatsORNULL", c("SolverStats", "NULL"))
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

# TODO: Get rid of this and just skip calling copy on Solution objects.
setMethod("copy", "Solution", function(object) {
  return(Solution(object@status, object@opt_val, object@primal_vars, object@dual_vars, object@attr))
})

#' @param x A \linkS4class{Solution} object.
#' @rdname Solution-class
setMethod("as.character", "Solution", function(x) {
  paste("Solution(status = ", x@status,
              ", opt_val = ", x@opt_val,
              ", primal_vars = (", paste(x@primal_vars, collapse = ", "),
              "), dual_vars = (", paste(x@dual_vars, collapse = ", "),
              "), attr = (", paste(x@attr, collapse = ", "), "))", sep = "")
})

setMethod("show", "Solution", function(object) {
  cat("Solution(", object@status, ", ",
                   object@opt_val, ", (",
                   paste(object@primal_vars, collapse = ", "), "), (",
                   paste(object@dual_vars, collapse = ", "), "), (",
                   paste(object@attr, collapse = ", "), "))", sep = "")
})

INF_OR_UNB_MESSAGE <- paste("The problem is either infeasible or unbounded, but the solver cannot tell which.",
                           "Disable any solver-specific presolve methods and re-solve to determine the precise problem status.",
                           "For GUROBI and CPLEX you can automatically perform this re-solve with the keyword argument solve(prob, reoptimize = TRUE, ...).")

# Factory function for infeasible or unbounded solutions.
failure_solution <- function(status, attr = NULL) {
  if(status %in% c(INFEASIBLE, INFEASIBLE_INACCURATE))
    opt_val <- Inf
  else if(status %in% c(UNBOUNDED, UNBOUNDED_INACCURATE))
    opt_val <- -Inf
  else
    opt_val <- NA_real_

  if(is.null(attr))
    attr <- list()
  if(status == INFEASIBLE_OR_UNBOUNDED)
    attr$message <- INF_OR_UNB_MESSAGE

  return(Solution(status, opt_val, list(), list(), attr))
}

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
#' @slot solution (Internal) Used internally to hold the problem solution
#' @slot cache (Internal) Used internally to hold cached matrix data.
#' @slot solver_cache (Internal) Used internally to hold cached solver data.
#' @slot size_metrics (Internal) Used internally to hold size metrics.
#' @slot solver_stats (Internal) Used internally to hold solver statistics.
#' @slot compilation_time (Internal) Used internally to hold the compilation time.
#' @slot solve_time (Internal) Used internally to hold the solve time.
#' @slot args (Internal) Used internally to hold the arguments (objective and constraints).
#' @name Problem-class
#' @aliases Problem
#' @rdname Problem-class
.Problem <- setClass("Problem", representation(objective = "Objective", constraints = "list", variables = "list", value = "numeric", status = "character", solution = "ANY", cache = "Cache", solver_cache = "list", size_metrics = "SizeMetricsORNULL", solver_stats = "SolverStatsORNULL", compilation_time = "numeric", solve_time = "numeric", args = "list"),
                    prototype(constraints = list(), value = NA_real_, status = NA_character_, solution = NULL, cache = Cache(), solver_cache = list(), size_metrics = NULL, solver_stats = NULL, compilation_time = NA_real_, solve_time = NA_real_, args = list()), contains = "Canonical")

#' @param objective A \linkS4class{Minimize} or \linkS4class{Maximize} object representing the optimization objective.
#' @param constraints (Optional) A list of \linkS4class{Constraint} objects representing constraints on the optimization variables.
#' @rdname Problem-class
#' @export
Problem <- function(objective, constraints = list()) {
  .Problem(objective = objective, constraints = constraints)
}

Problem.validate_constraint <- function(constraint) {
  if(is(constraint, "Constraint"))
    return(constraint)
  else if(is(constraint, "logical")) {
    # Replace TRUE or FALSE values with equivalent Expressions.
    if(constraint)
      return(Constant(0) <= Constant(1))
    else
      return(Constant(1) <= Constant(0))
  } else
    stop("Problem has an invalid constraint of type ", class(constraint))
}

setMethod("initialize", "Problem", function(.Object, ..., objective, constraints = list(), value = NA_real_, status = NA_character_, solution = NULL, cache = Cache(), solver_cache = list(), size_metrics = NULL, solver_stats = NULL, compilation_time = NA_real_, solve_time = NA_real_, args = list()) {
  if(is.null(constraints))
    constraints <- list()

  # Check that objective is Minimize or Maximize.
  if(!is(objective, "Minimize") && !is(objective, "Maximize"))
    stop("Problem objective must be Minimize or Maximize")

  # Constraints and objective are immutable
  .Object@objective <- objective
  # Raise warning if objective has too many subexpressions.
  if(node_count(.Object@objective) >= MAX_NODES)
    warning("Objective constraints too many subexpressions. Consider vectorizing your CVXR code to speed up compilation")

  .Object@constraints <- lapply(constraints, Problem.validate_constraint)
  # Raise warning if constraint has too many subexpressions.
  for(i in seq_along(.Object@constraints)) {
    constraint <- .Object@constraints[[i]]
    if(node_count(constraints) >= MAX_NODES)
      warning("Constraint ", i, " contains too many subexpressions. Consider vectorizing your CVXR code to speed up compilation")
  }

  .Object@value <- NA_real_
  .Object@status <- NA_character_
  .Object@solution <- NULL
  .Object@cache <- Cache()
  .Object@solver_cache <- list()

  # Information about the dimensions of the problem and its constituent parts.
  .Object@size_metrics <- NULL

  # Benchmarks reported by the solver.
  .Object@solver_stats <- NULL
  .Object@compilation_time <- NA_real_
  .Object@solve_time <- NA_real_
  .Object@args <- list(.Object@objective, .Object@constraints)
  .Object
})

#' @param object A \linkS4class{Problem} object.
#' @describeIn Problem The value from the last time the problem was solved (or NA if not solved).
setMethod("value", "Problem", function(object) {
  if(is.na(object@value))
    return(NA_real_)
  else
    return(intf_scalar_value(object@value))
})

#' @describeIn Problem The status from the last time the problem was solved; one of optimal, infeasible, or unbounded (with or without suffix inaccurate).
setMethod("status", "Problem", function(object) { object@status })

#' @describeIn Problem The solution from the last time the problem was solved.
setMethod("solution", "Problem", function(object) { object@solution })

#' @describeIn Problem The objective of the problem. Note that the objective cannot be reassigned after creation, and modifying the objective after creation will result in undefined behavior.
setMethod("objective", "Problem", function(object) { object@objective })

#' @describeIn Problem A list of the constraints of the problem. Note that constraints cannot be reassigned, appended to, or otherwise modified after creation, except through parameters.
setMethod("constraints", "Problem", function(object) { object@constraints })

#' @describeIn Problem Expose all parameters as a named list.
setMethod("param_list", "Problem", function(object) {
  res <- list()
  for(parm in parameters(object))
    res[[name(parm)]] <- parm
  return(res)
})

#' @describeIn Problem Expose all variables as a named list.
setMethod("var_list", "Problem", function(object) {
  res <- list()
  for(var in variables(problem))
    res[[name(var)]] <- var
  return(res)
})

#' @param dpp A logical value indicating whether to enforce the disciplined parametrized programming (DPP) ruleset; only relevant when the problem involves Parameters. DPP is a mild restriction of DCP. When a problem involving Parameters is DPP, subsequent solves can be much faster than the first one. For more information, consult the documentation at https://www.cvxpy.org/tutorial/advanced/index.html#disciplined-parametrized-programming
#' @describeIn Problem A logical value indicating whether the problem satisfies DCP rules.
#' @examples
#' x <- Variable(2)
#' p <- Problem(Minimize(p_norm(x, 2)), list(x >= 0))
#' is_dcp(p)
setMethod("is_dcp", "Problem", function(object, dpp = FALSE) {
  all(sapply(c(object@constraints, list(object@objective)), function(expr) { is_dcp(expr, dpp) }))
})

#' @describeIn Problem A logical value indicating whether the problem statisfies DGP rules.
setMethod("is_dgp", "Problem", function(object, dpp = FALSE) {
  all(sapply(c(object@constraints, list(object@objective)), function(expr) { is_dgp(expr, dpp) }))
})

#' @describeIn Problem Does the problem satisfy the DQCP rules?
setMethod("is_dqcp", "Problem", function(object) {
  all(sapply(c(object@constraints, list(object@objective)), is_dqcp))
})

#' @param context A string indicating whether to check DPP-compliance for DCP (\code{context = "dcp"}) or DGP (\code{context = "dgp"}).
#' @describeIn Problem Does the problem satisfy DPP rules?
setMethod("is_dpp", "Problem", function(object, context = "dcp") {
  if(tolower(context) == "dcp")
    return(is_dcp(object, dpp = TRUE))
  else if(tolower(context) == "dgp")
    return(is_dgp(object, dpp = TRUE))
  else
    stop("Unsupported context ", context)
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

#' @describeIn Problem logical value indicating whether the problem is a mixed integer program.
setMethod("is_mixed_integer", "Problem", function(object) {
  any(sapply(variables(object), function(v) { v@attributes$boolean || v@attributes$integer }))
})

#' @describeIn Problem List of \linkS4class{Variable} objects in the problem.
setMethod("variables", "Problem", function(object) {
  vars_ <- variables(object@objective)
  for(constr in object@constraints)
    vars_ <- c(vars_, variables(constr))
  unique_list(vars_)
})

#' @describeIn Problem List of \linkS4class{Parameter} objects in the problem.
setMethod("parameters", "Problem", function(object) {
  params <- parameters(object@objective)
  for(constr in object@constraints)
    params <- c(params, parameters(constr))
  unique_list(params)
})

#' @describeIn Problem List of \linkS4class{Constant} objects in the problem.
setMethod("constants", "Problem", function(object) {
  const_list <- list()
  constants_ <- constraints(object@objective)
  for(constr in object@constraints)
    constants_ <- c(constants_, constants(constr))
  # for(constant in constants_)
  #   const_dict[[as.character(id(constant))]] <- constant
  # unname(const_dict)
  unique_list(constants)
})

#' @describeIn Problem List of \linkS4class{Atom} objects in the problem.
setMethod("atoms", "Problem", function(object) {
  atoms_ <- atoms(object@objective)
  for(constr in object@constraints)
    atoms_ <- c(atoms_, atoms(constr))
  unique_list(atoms_)
})

#' @describeIn Problem Information about the size of the problem.
setMethod("size_metrics", "Problem", function(object) {
  if(is.null(object@size_metrics))
    object@size_metrics <- SizeMetrics(object)
  list(object, object@size_metrics)   # Return object because slot may have been modified.
})

#' @describeIn Problem Additional information returned by the solver.
setMethod("solver_stats", "Problem", function(object) { object@solver_stats })

# Compiles and solves the problem using the specified method.
# solve.Problem <- function(object, method, ...) {
#  if(missing(method))
#    .solve(object, ...)
#  else {
#    func <- Problem.REGISTERED_SOLVE_METHODS[[method]]
#    func(object, ...)
#  }
# }

# Adds a solve method to the Problem class.
# Problem.register_solve <- function(name, func) {
#  Problem.REGISTERED_SOLVE_METHODS[[name]] <- func
# }

#'
#' @param object A \linkS4class{Problem} class.
#' @param solver A string indicating the solver that the problem data is for. Call \code{installed_solvers()} to see all available.
#' @param gp A logical value indicating whether to parse the problem as a disciplined geometric program (DGP) instead of a disciplined convex program (DCP).
#' @param enforce_dpp If TRUE, a DPPError will be thrown when trying to parse a non-DPP problem (instead of just a warning). Defaults to FALSE.
#' @param ignore_dpp If TRUE, DPP problems will be treated as non-DPP, which may speed up compilation. Defaults to FALSE.
#' @param canon_backend Specifies which backend to use for canonicalization, which can affect compilation time. Defaults to 'CPP'.
#' @param verbose If TRUE, print verbose output related to problem compilation.
#' @param solver_opts A list of options that will be passed to the specific solver. In general, these options will override any default settings imposed by CVXR.
#' @describeIn Problem Get the problem data passed to the specified solver.
setMethod("get_problem_data", "Problem", function(object, solver, gp = FALSE, enforce_dpp = FALSE, ignore_dpp = FALSE, verbose = FALSE, canon_backend = NA_character_, solver_opts = NULL) {
  # Invalid DPP setting.
  # Must be checked here to avoid cache issues.
  if(enforce_dpp && ignore_dpp)
    stop("DPPError: Cannot set enforce_dpp = TRUE and ignore_dpp = TRUE)")

  start_time <- Sys.time()

  # Cache includes ignore_dpp because it alters compilation.
  key <- make_key(object@cache, solver, gp, ignore_dpp)
  if(length(key) != length(object@cache@key) || !all(mapply(function(x, y) { x == y }, key, object@cache@key))) {
    object@cache <- invalidate(object@cache)
    solving_chain <- .construct_chain(object, solver = solver, gp = gp, enforce_dpp = enforce_dpp, ignore_dpp = ignore_dpp, canon_backend = canon_backend, solver_opts = solver_opts)
    object@cache@key <- key
    object@cache@solving_chain <- solving_chain
    object@solver_cache <- list()
  } else
    solving_chain <- object@cache@solving_chain

  if(verbose) {
    divider <- paste(rep("-", 79), collapse = "")
    cat(paste(divider, sprintf("%45.11s", "Compilation"), divider, sep = "\n"))
  }

  if(!is.null(object@cache@param_prog)) {
    # Fast path, bypasses application of reductions.
    if(verbose)
      print("Using cached ASA map, for faster compilation (bypassing reduction chain)")

    if(gp) {
      dgp2dcp <- get(object@cache@solving_chain, "Dgp2Dcp")
      # Parameters in the param cone prog are the logs of parameters in the original problem (with one exception: parameters appearning as exponents (in Power and GmatMul atoms) are unchanged).
      old_params_to_new_params <- dgp2dcp@canon_methods@parameters
      for(param in parameters(object)) {
        pid <- as.character(id(param))
        if(pid %in% names(old_params_to_new_params))
          value(old_params_to_new_params[[pid]]) <- log(value(param))   # TODO: Need to propagate new parameter values to param_prog.
      }
    }

    tmp <- perform(solving_chain@solver, object@cache@param_prog)
    data <- tmp[[1]]
    solver_inverse_data <- tmp[[2]]
    inverse_data <- c(object@cache@inverse_data, list(solver_inverse_data))
    object@compilation_time <- Sys.time() - start_time
    object@compilation_time <- as.numeric(object@compilation_time)

    if(verbose)
      print(paste("Finished problem compilation (took", sprintf("%.3e", object@compilation_time), "seconds"))
  } else {
    if(verbose) {
      solver_name <- name(solving_chain@reductions[[length(solving_chain@reductions)]])
      reduction_chain_str <- paste(sapply(solving_chain@reductions, class), collapse = " -> ")
      print(paste("Compiling problem (target solver = ", solver_name, ")", sep = ""))
      print(paste("Reduction chain:", reduction_chain_str))
    }

    tmp <- perform(solving_chain, object, verbose)
    data <- tmp[[1]]
    inverse_data <- tmp[[2]]
    save_to_cache <- (is(data, "list") && PARAM_PROB %in% names(data) && !any(sapply(solving_chain@reductions, function(reduction) { is(reduction, "EvalParams") })))
    object@compilation_time <- Sys.time() - start_time
    object@compilation_time <- as.numeric(object@compilation_time)

    if(verbose)
      print(paste("Finished problem compilation (took", sprintf("%.3e", object@compilation_time), "seconds"))

    if(safe_to_cache) {
      if(verbose && !is.null(parameters(object)) && length(parameters(object)) > 0)
        print("(Subsequent compilations of this problem, using the same arguments, should take less time)")
      object@cache@param_prog <- data[[PARAM_PROB]]
      # The last datum in inverse_data corresponds to the solver, so we shouldn't cache it.
      if(length(inverse_data) <= 1)
        object@cache@inverse_data <- NULL
      else
        object@cache@inverse_data <- inverse_data[seq(length(inverse_data) - 1)]
    }
  }

  return(list(data = data, solving_chain = solving_chain, inverse_data = inverse_data))
})

# Find candidate solvers for the current problem. If solver is not NA, it checks if the specified solver is compatible with the problem passed.
.find_candidate_solvers <- function(object, solver = NA_character_, gp = FALSE) {
  candidates <- list(qp_solvers = list(), conic_solvers = list())
  if(is(solver, "Solver"))
    return(.add_custom_solver_candidates(solver))

  if(!is.na(solver)) {
    if(!(solver %in% INSTALLED_SOLVERS))
      stop("The solver ", solver, " is not installed")
    if(solver %in% CONIC_SOLVERS)
      candidates$conic_solvers <- c(candidates$conic_solvers, list(solver))
    if(solver %in% QP_SOLVERS)
      candidates$qp_solvers <- c(candidates$qp_solvers, list(solver))
  } else {
    candidates$qp_solvers <- INSTALLED_SOLVERS[INSTALLED_SOLVERS %in% QP_SOLVERS]
    candidates$conic_solvers <- list()
    # ECOS_BB can only be called explicitly.
    for(slv in INSTALLED_SOLVERS) {
      if(slv %in% CONIC_SOLVERS && slv != ECOS_BB)
        candidates$conic_solvers <- c(candidates$conic_solvers, list(slv))
    }
  }

  # If gp, we must have only conic solvers.
  if(gp) {
    if(!is.na(solver) && !(solver %in% CONIC_SOLVERS))
      stop("When gp = TRUE, solver must be a conic solver (received ", solver, "); try calling solve() with solver = 'ECOS'")
    else if(is.na(solver))
      candidates$qp_solvers <- list()   # No QP solvers allowed.
  }

  if(is_mixed_integer(object)) {
    # ECOS_BB must be called explicitly.
    if(identical(INSTALLED_MI_SOLVERS, list(ECOS_BB)) && solver != ECOS_BB)
      stop("You need a mixed-integer solver for this model. Refer to the documentation
            https://www.cvxpy.org/tutorial/advanced/index.html#mixed-integer-programs
            for discussion on this topic.

            Quick fix 1: if you install the python package CVXOPT (pip install cvxopt),
            then CVXR can use the open-source mixed-integer linear programming
            solver `GLPK`. If your problem is nonlinear then you can install SCIP
            (pip install pyscipopt).

            Quick fix 2: you can explicitly specify solver = 'ECOS_BB'. This may result
            in incorrect solutions and is not recommended.")

    # TODO: Provide a useful error message when the problem is nonlinear, but the only
    # installed mixed-integer solvers are MILP solvers (e.g., GLPK_MI).
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

# Returns a list of candidate solvers where custom_solver is the only potential option.
.add_custom_solver_candidates <- function(object, custom_solver) {
  if(name(custom_solver) %in% SOLVERS)
    stop("Custom solvers must have a different name than the officially supported ones")

  candidates <- list("qp_solvers" = list(), "conic_solvers" = list())
  if(!is_mixed_integer(object) || mip_capable(custom_solver)) {
    if(is(custom_solver, "QpSolver")) {
      SOLVER_MAP_QP[[name(custom_solver)]] <- custom_solver
      candidates$qp_solvers <- list(name(custom_solver))
    } else if(is(custom_solver, "ConicSolver")) {
      SOLVER_MAP_CONIC[[name(custom_solver)]] <- custom_solver
      candidates$conic_solvers <- list(name(custom_solver))
    }
  }
  return(candidates)
}

# Construct the chains required to reformulate and solve the problem.
# In particular, this function 1) finds the candidates solvers, 2) constructs the solving chain that performs the numeric reductions and solves the problem.
.construct_chain <- function(object, solver = NA_character_, gp = FALSE, enforce_dpp = FALSE, ignore_dpp = FALSE, canon_backend = NA_character_, solver_opts = NULL) {
  candidate_solvers <- .find_candidate_solvers(solver = solver, gp = gp)
  .sort_candidate_solvers(object, candidate_solvers)
  return(construct_solving_chain(object, candidate_solvers, gp = gp, enforce_dpp = enforce_dpp, ignore_dpp = ignore_dpp, canon_backend = canon_backend, solver_opts = solver_opts))
}

# Sorts candidate solvers lists according to CONIC_SOLVERS/QP_SOLVERS.
Problem.sort_candidate_solvers <- function(solvers) {
  if(length(solvers$conic_solvers) > 1) {
    idx <- match(solvers$conic_solvers, CONIC_SOLVERS)
    solvers$conic_solvers <- CONIC_SOLVERS[sort(idx, decreasing = FALSE)]
  }
  if(length(solvers$qp_solvers) > 1) {
   idx <- match(solvers$qp_solvers, QP_SOLVERS)
   solvers$qp_solvers <- QP_SOLVERS[sort(idx, decreasing = FALSE)]
  }
  return(solvers)
}

.invalidate_cache <- function(object) {
  object@cache_key <- NA_character_
  object@solving_chain <- NULL
  object@param_prog <- NULL
  object@inverse_data <- NULL
  return(object)
}

#' @docType methods
#' @rdname Problem This internal method solves a DCP compliant optimization problem. It saves the values of primal and dual variables in the Variable and Constraint objects, respectively.
#' @export
setMethod("psolve", "Problem", function(object, solver = NA_character_, ignore_dcp = FALSE, warm_start = TRUE, verbose = FALSE, gp = FALSE, qcp = FALSE,
                                        requires_grad = FALSE, enforce_dpp = FALSE, ignore_dpp = FALSE, canon_backend = NA_character_,
                                        parallel = FALSE, feastol = NULL, reltol = NULL, abstol = NULL, num_iter = NULL, low = NA_real_, high = NA_real_, ...) {
  # TODO: Need to update this function.
  # if(parallel)
  #   stop("Unimplemented")

  if(verbose) {
    divider <- paste(rep("=", 79), collapse = "")
    version <- packageVersion("CVXR")
    cat(paste(divider, sprintf("%38.4s", "CVXR"), sprintf("%39.6s", paste("v", version, sep = "")), divider, sep = "\n"))
  }

  for(parameter in parameters(object)) {
    if(is.na(value(parameter)))
      stop("A Parameter (whose name is ", name(parameter), " ) does not have a value associated with it; all Parameter objects must have values before solving a problem")
  }

  if(verbose) {
    n_variables <- sum(sapply(variables(object), function(v) { prod(dim(v)) }))
    n_parameters <- sum(sapply(parameters(object), function(p) { prod(dim(p)) }))
    print(paste("Your problem has", n_variables, "variables,", length(object@constraints), "constraints, and", n_parameters, "parameters"))

    curvatures <- c()
    if(is_dcp(object))
      curvatures <- c(curvatures, "DCP")
    if(is_dgp(object))
      curvatures <- c(curvatures, "DGP")
    if(is_dqcp(object))
      curvatures <- c(curvatures, "DQCP")
    print(paste("It is compliant with the following grammars:", paste(curvatures, collapse = ", ")))

    if(n_parameters == 0)
      print("(If you need to solve this problem multiple times, but with different data, consider using parameters)")
    print("CVXR will first compile your problem;, then, it will invoke a numerical solver to obtain a solution")
  }

  if(requires_grad) {
    dpp_context <- ifelse(gp, "dgp", "dcp")
    if(qcp)
      stop("Cannot compute gradients of DQCP problems")
    else if(!is_dpp(object, dpp_context))
      stop("Problem is not DPP (when requires_grad = TRUE, problem must be DPP)")
    else if(!is.na(solver) && !is(solver, "SCS") && !is(solver, "DIFFCP"))
      stop("When requires_grad = TRUE, the only supported solver is SCS")
    else if(!("DIFFCP" %in% INSTALLED_SOLVERS))
      stop("The R library diffcp must be installed to differentiate through problems")
    else
      solver <- DIFFCP
  } else {
    if(gp && qcp)
      stop("At most one of gp and qcp can be TRUE")
    if(qcp && !is_dcp(object)) {
      if(!is_dqcp(object))
        stop("The problem is not DQCP")
      if(verbose)
        stop("Reducing DQCP problem to a one-parameter family of DCP problems, for bisection")
      reductions <- list(Dqcp2Dcp())
      start_time <- Sys.time()
      if(inherits(object@objective, "Maximize")) {
        reductions <- c(list(FlipObjective), reductions)
        # FlipObjective flips the sign of the objective
        low_in <- low
        high_in <- high
        if(!is.na(high_in))
          low <- high_in*-1
        if(!is.na(low_in))
          high <- low_in*-1
      }

      chain <- Chain(problem = object, reductions = reductions)
      soln <- bisect(reduce(chain), solver = solver, verbose = verbose, low = low, high = high, verbose = verbose, ...)
      object <- unpack_problem(object, retrieve(chain, soln))
      return(object@value)
    }
  }

  kwargs <- c(list(feastol = feastol, reltol = reltol, abstol = abstol, num_iter = num_iter, low = low, high = high), list(...))
  tmp <- get_problem_data(object, solver, gp, enforce_dpp, ignore_dpp, verbose, canon_backend, kwargs)
  data <- tmp[[1]]
  solving_chain <- tmp[[2]]
  inverse_data <- tmp[[3]]

  if(verbose) {
    divider <- paste(rep("-", 79), collapse = "")
    cat(paste(divider, sprintf("%45.16s", "Numerical solver"), divider, sep = "\n"))
    solver_name <- name(solving_chain@reductions[[length(solving_chain@reductions)]])
    print(paste("Invoking solver", solver_name, "to obtain a solution"))
  }

  start_time <- Sys.time()
  solution <- solve_via_data(solving_chain, object, data, warm_start, verbose, kwargs)
  end_time <- Sys.time()
  object@solve_time <- as.numeric(end_time - start_time)
  object <- unpack_results(object, solution, solving_chain, inverse_data)

  if(verbose) {
    divider <- paste(rep("-", 79), collapse = "")
    cat(paste(divider, sprintf("%40.7s", "Summary"), divider, sep = "\n"))
    print(paste("Problem status:", object@status))
    print(paste("Optimal value:", object@value))
    print(paste("Compilation took", object@compilation_time, "seconds"))
    print(paste("Solver (including time spent in interface) took", object@solve_time, "seconds"))
  }
  return(list(object, object@value))
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

#' @describeIn Problem Compute the gradient of a solution with respect to parameters.
setMethod("backward", "Problem", function(object) {
  if(!(DIFFCP %in% object@solver_cache))
    stop("backward can only be called after calling solve with requires_grad = TRUE")
  else if(!(object@status %in% SOLUTION_PRESENT))
    stop("Backpropagating through infeasible/unbounded problems is not yet supported. Please file an issue on Github if you need this feature")

  # TODO: Backpropagate through dual variables as well.
  backward_cache <- object@solver_cache[[DIFFCP]]
  DT <- backward_cache$DT
  zeros <- matrix(0, nrow = nrow(backward_cache$s), ncol = ncol(backward_cache$s))
  del_vars <- list()

  gp <- gp(object@cache)
  for(variable in variables(object)) {
    vid_char <- as.character(id(variable))
    if(is.na(variable@gradient))
      del_vars[[vid_char]] <- matrix(1, nrow = nrow(variable), ncol = ncol(variable))
    else
      del_vars[[vid_char]] <- as.vector(variable@gradient)

    if(gp) {
      # x_gp <- exp(x_cone_program),
      # dx_gp/d x_cone_program <- exp(x_cone_program) <- x_gp
      del_vars[[vid_char]] <- del_vars[[vid_char]]*value(variable)
    }
  }

  dx <- split_adjoint(object@cache@param_prog, del_vars)
  start_time <- Sys.time()
  tmp <- DT(dx, zeros, zeros)
  dA <- tmp[[1]]
  db <- tmp[[2]]
  dc <- tmp[[3]]
  end_time <- Sys.time()
  backward_cache$DT_TIME <- end_time - start_time
  backward_cache$DT_TIME <- as.numeric(backward_cache$DT_TIME)
  dparams <- apply_param_jac(object@cache@param_prog, dc, -dA, db)

  if(!gp) {
    for(i in seq_along(object@parameters)) {
      pid_char <- as.character(id(object@parameters[[i]]))
      object@parameters[[i]]@gradient <- dparams[[pid_char]]
    }
  } else {
    dgp2dcp <- get(object@cache@solving_chain, "Dgp2Dcp")
    old_params_to_new_params <- dgp2dcp@canon_methods@parameters
    for(i in seq_along(object@parameters)) {
      param <- object@parameters[[i]]
      # Note: If param is an exponent in a Power or GmatMul atom, then the
      # parameter passes through unchanged to the DCP program; if the param is
      # also used elsewhere (not as an exponent), then param will also be in
      # ol_params_to_new_params. Therefore,
      # param@gradient = dparams[[as.character(id(param))]] (or 0) + 1/param*dparams[[as.character(id(new_param))]]
      #
      # Note that id(param) is in dparams if and only if param was used as an
      # exponent (because this means that the parameter entered the DCP problem
      # unchanged).
      pid_char <- as.character(id(param))
      if(pid_char %in% names(dparams))
        grad <- dparams[[pid_char]]
      else
        grad <- 0.0

      if(pid_char %in% names(old_params_to_new_params)) {
        new_param <- old_params_to_new_params[[pid_char]]
        # value(new_param) == log(param), apply chain rule.
        grad <- grad + (1.0 / value(param)) * dparams[[id(new_param)]]
      }
      object@parameters[[i]]@gradient <- grad
    }
  }
  return(object)
})

#' @describeIn Problem Apply the derivative of the solution map to perturbations in the parameters.
setMethod("derivative", "Problem", function(object) {
  if(!(DIFFCP %in% object@solver_cache))
    stop("derivative can only be called after calling solve with requires_grad = TRUE")
  else if(!(object@status %in% SOLUTION_PRESENT))
    stop("Differentiating through infeasible/unbounded problems is not yet supported. Please file an issue on Github if you need this feature")

  # TODO: Forward differentiate dual variables as well.
  backward_cache <- object@solver_cache[[DIFFCP]]
  param_prog <- object@cache@param_prog
  D <- backward_cache$D
  param_deltas <- list()

  gp <- gp(object@cache)
  if(gp)
    dgp2dcp <- get(object@cache@solving_chain, "Dgp2Dcp")

  if(is.null(parameters(object)) || length(parameters(object)) == 0) {
    for(i in seq_along(object@variables)) {
      vdim <- dim(object@variables[[i]])
      object@variables[[i]]@delta <- matrix(0, nrow = vdim[1], ncol = vdim[2])
    }
    return(object)
  }

  for(param in parameters(object)) {
    pid <- id(param)
    pid_char <- as.character(pid)

    if(is.na(param@delta))
      delta <- matrix(0, nrow = nrow(param), ncol = ncol(param))
    else
      delta <- param@delta
    if(gp) {
      if(pid_char %in% names(dgp2dcp@canon_methods@parameters))
        new_param_id <- id(dgp2dcp@canon_methods@parameters[[pid_char]])
      else
        new_param_id <- pid
      param_deltas[[as.character(new_param_id)]] <- 1.0/value(param) * as.vector(delta)
      if(as.character(id(param)) %in% names(param_prog@param_id_to_col)) {
        # Here, param generated a new parameter and also passed through to the param cone prog unchanged (because it was an exponent of a power).
        param_deltas[[pid_char]] <- as.vector(delta)
      }
    } else
      param_deltas[[pid_char]] <- as.vector(delta)
  }

  tmp <- apply_parameters(param_prog, param_deltas, zero_offset = TRUE)
  dc <- tmp[[1]]
  dA <- tmp[[3]]
  db <- tmp[[4]]

  start_time <- Sys.time()
  dx <- D(-dA, db, dc)[[1]]
  end_time <- Sys.time()
  backward_cache$D_TIME <- end_time - start_time
  backward_cache$D_TIME <- as.numeric(backward_cache$D_TIME)
  dvars <- split_solution(param_prog, dx, sapply(variables(object), id))

  for(i in seq_along(object@variables)) {
    vid_char <- as.character(id(object@variables[[i]]))
    object@variables[[i]]@delta <- dvars[[vid_char]]
    if(gp) {
      # x_gp = exp(x_cone_program),
      # dx_gp/d x_cone_program = exp(x_cone_program) = x_gp
      var_val <- value(object@variables[[i]])
      object@variables[[i]]@delta <- object@variables[[i]]@delta*var_val
    }
  }
  return(object)
})

.clear_solution <- function(object) {
  for(i in seq_along(object@variables)) {
    # value(object@variables[[i]]) <- NA_real_
    object@variables[[i]] <- save_value(object@variables[[i]], NA_real_)
  }

  for(i in seq_along(object@constraints)) {
    for(j in seq_along(object@constraints[[i]]@dual_variables)) {
      # value(object@constraints[[i]]@dual_variables[[j]]) <- NA_real_
      object@constraints[[i]]@dual_variables[[j]] <- save_value(object@constraints[[i]]@dual_variables[[j]], NA_real_)
    }
  }

  object@value <- NA_real_
  object@status <- NA_character_
  object@solution <- NULL
  return(object)
}

setMethod("unpack_problem", signature(object = "Problem", solution = "Solution"), function(object, solution) {
  # if(solution@status %in% SOLUTION_PRESENT) {
  #   for(i in seq_along(object@variables)) {
  #     vid_char <- as.character(id(object@variables[[i]]))
  #     object@variables[[i]] <- save_value(object@variables[[i]], solution@primal_vars[[vid_char]])
  #   }
  #
  #   for(i in seq_along(object@constraints)) {
  #     cid_char <- as.character(id(object@constraints[[i]]))
  #     if(cid_char %in% names(solution@dual_vars))
  #       object@constraints[[i]] <- save_dual_value(object@constraints[[i]], solution@dual_vars[[cid_char]])
  #   }
  #
  #   # Eliminate confusion of problem@value versus objective@value
  #   object@value <- object@objective@value
  # } else if(solution@status %in% INF_OR_UNB) {
  #   for(i in seq_along(object@variables))
  #     object@variables[[i]] <- save_value(object@variables[[i]], NA_real_)
  #
  #   for(i in seq_along(object@constraints)) {
  #     for(j in seq_along(object@constraints[[i]]@dual_variables))
  #       object@constraints[[i]]@dual_variables[[j]] <- save_value(object@constraints[[i]]@dual_variables[[j]], NA_real_)
  #   }
  #   object@value <- solution@opt_val
  # } else
  #   stop("Cannot unpack invalid solution")
  #
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
        if(is.null(dim(val)) || all(dim(val) == 1))
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
#' @param solution A \linkS4class{Solution} object.
#' @param chain The corresponding solving \linkS4class{Chain}.
#' @param inverse_data A \linkS4class{InverseData} object or list containing data necessary for the inversion.
#' @docType methods
setMethod("unpack_results", "Problem", function(object, solution, chain, inverse_data) {
  solution <- invert(chain, solution, inverse_data)

  if(solution@status %in% INACCURATE)
    warning("Solution may be inaccurate. Try another solver, adjusting the solver settings, or solve with verbose = TRUE for more information")
  if(solution@status == INFEASIBLE_OR_UNBOUNDED)
    warning(paste("The problem is either infeasible or unbounded, but the solver cannot tell which. Disable any solver-specific presolve methods and re-solve to determine the precise problem status.",
                  "For GUROBI and CPLEX, you can automatically perform this re-solve with the keyword argument solve(prob, reoptimize = TRUE, ...)", sep = "\n"))
  if(solution@status %in% ERROR)
    stop("Solver failed. Try another solver, or solve with verbose = TRUE for more information")

  object <- unpack_problem(object, solution)
  object@solver_stats <- SolverStats(object@solution@attr, name(chain@solver))
  return(object)

  # results <- unpack_problem(object, solution)
  # solver_stats <- SolverStats(solution@attr, name(chain@solver))
  # return(c(results, solver_stats))
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
        if(inherits(constr, constr_types))
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

        ## Convert matrix back to scalar/vector when appropriate.
        # if(var@.is_vector)
        #  value <- as.vector(value)
        value
    })
    names(result) = sapply(variables, function(x) as.character(id(x)))
    result
}

setMethod("as.character", "Problem", function(x) {
  if(length(object@constraints) == 0)
    return(as.character(object@objective))
  else {
    subect_to <- "subject to "
    lines <- c(as.character(object@objective), paste(subject_to, object@constraints[[1]], sep = ""))
    if(length(object@constraints) > 1) {
      for(i in seq(2, length(object@constraints))) {
        constr <- object@constraints[[i]]
        lines <- c(lines, paste(rep(" ", length(subject_to)), as.character(constr), collapse = ""))
      }
    }
    return(paste(lines,  collapse = "\n"))
  }
})

setMethod("show", "Problem", function(object) {
  cat("Problem(", as.character(object@objective), ", (", paste(sapply(object@constraints, as.character), ", "), "))", sep = "")
})

setMethod("is_constant", "Problem", function(object) { FALSE })

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
  Problem(objective = e1@objective + e2@objective, constraints = unique_list(c(e1@constraints, e2@constraints)))
})

#' @rdname Problem-arith
setMethod("-", signature(e1 = "Problem", e2 = "numeric"), function(e1, e2) { e1 + (-e2) })

#' @rdname Problem-arith
setMethod("-", signature(e1 = "numeric", e2 = "Problem"), function(e1, e2) { if(length(e1) == 1 && e1 == 0) -e2 else stop("Unimplemented") })

#' @rdname Problem-arith
setMethod("-", signature(e1 = "Problem", e2 = "Problem"), function(e1, e2) {
  Problem(objective = e1@objective - e2@objective, constraints = unique_list(c(e1@constraints, e2@constraints)))
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

#'
#' The ParamProb class.
#'
#' This virtual class represents parametrized problems.
#'
#' Parameterized problems are produced during the first canonicalization
#' and allow canonicalization to be short-circuited for future solves.
#'
#' @name ParamProb-class
#' @aliases ParamProb
#' @rdname ParamProb-class
ParamProb <- setClass("ParamProb", contains = "VIRTUAL")

#' @param object A \linkS4class{ParamProb} object.
#' @describeIn ParamProb Is the problem mixed-integer?
setMethod("is_mixed_integer", "ParamProb", function(object) { stop("Unimplemented") })

#' @param id_to_param_value (Optional) List mapping parameter IDs to values.
#' @param zero_offset (Optional) If TRUE, zero out the constant offset in the parameter vector.
#' @param keep_zeros (Optional) If TRUE, store explicit zeros in A where parameters are affected.
#' @describeIn ParamProb Returns A, b after applying parameters (and reshaping).
setMethod("apply_parameters", "ParamProb", function(object, id_to_param_value = NULL, zero_offset = FALSE, keep_zeros = FALSE) {
  stop("Unimplemented")
})

#'
#' A Partial Optimization Problem
#'
#' @slot opt_vars The variables to optimize over.
#' @slot dont_opt_vars The variables to not optimize over.
#' @name PartialProblem-class
#' @aliases PartialProblem
#' @rdname PartialProblem-class
.PartialProblem <- setClass("PartialProblem", representation(.prob = "Problem", opt_vars = "list", dont_opt_vars = "list", solver = "ANY", .solve_kwargs = "list"),
                                              prototype(opt_vars = list(), dont_opt_vars = list(), solver = NA_character_, .solve_kwargs = list()),
                             contains = "Expression")

PartialProblem <- function(prob, opt_vars, dont_opt_vars, solver, ...) {
  .PartialProblem(.prob = prob, opt_vars = opt_vars, dont_opt_vars = dont_opt_vars, solver = solver, .solve_kwargs = list(...))
}

setMethod("initialize", "PartialProblem", function(.Object, ..., .prob, opt_vars = list(), dont_opt_vars = list(), solver = NA_character_, .solve_kwargs = list()) {
  .Object@.prob <- .prob
  .Object@opt_vars <- opt_vars
  .Object@dont_opt_vars <- dont_opt_vars
  .Object@solver <- solver
  .Object@.solve_kwargs <- .solve_kwargs
  callNextMethod(.Object, ..., args = list(.prob))
})

#' @describeIn PartialProblem Returns info needed to reconstruct the expression besides the args.
setMethod("get_data", "PartialProblem", function(object) { list(object@opt_vars, object@dont_opt_vars, object@solver) })

#' @describeIn PartialProblem Is the expression constant?
setMethod("is_constant", "PartialProblem", function(object) {
  length(variables(object@args[[1]])) == 0
})

#' @describeIn PartialProblem Is the expression convex?
setMethod("is_convex", "PartialProblem", function(object) {
  is_dcp(object@args[[1]]) && inherits(object@args[[1]]@objective, "Minimize")
})

#' @describeIn PartialProblem Is the expression concave?
setMethod("is_concave", "PartialProblem", function(object) {
  is_dcp(object@args[[1]]) && inherits(object@args[[1]]@objective, "Maximize")
})

#' @describeIn PartialProblem Is the expression a disciplined parameterized expression?
setMethod("is_dpp", "PartialProblem", function(object, context = "dcp") {
  if(tolower(context) %in% c("dcp", "dgp"))
    is_dpp(object@args[[1]], context)
  else
    stop("Unsupported context ", context)
})

#' @describeIn PartialProblem Is the expression log-log convex?
setMethod("is_log_log_convex", "PartialProblem", function(object) {
  is_dgp(object@args[[1]]) && inherits(object@args[[1]]@objective, "Minimize")
})

#' @describeIn PartialProblem Is the expression log-log concave?
setMethod("is_log_log_concave", "PartialProblem", function(object) {
  is_dgp(object@args[[1]]) && inherits(object@args[[1]]@objective, "Maximize")
})

#' @describeIn PartialProblem Is the expression nonnegative?
setMethod("is_nonneg", "PartialProblem", function(object) {
  is_nonneg(object@args[[1]]@objective@args[[1]])
})

#' @describeIn PartialProblem Is the expression nonpositive?
setMethod("is_nonpos", "PartialProblem", function(object) {
  is_nonpos(object@args[[1]]@objective@args[[1]])
})

#' @describeIn PartialProblem Is the expression imaginary?
setMethod("is_imag", "PartialProblem", function(object) { FALSE })

#' @describeIn PartialProblem Is the expression complex valued?
setMethod("is_complex", "PartialProblem", function(object) { FALSE })

#' @describeIn PartialProblem Returns the (row, col) dimensions of the expression.
setMethod("dim", "PartialProblem", function(x) { c(1,1) })

setMethod("name", "PartialProblem", function(x) { cat("PartialProblem(", x@args[[1]], ")") })

#' @describeIn PartialProblem Returns the variables in the problem.
setMethod("variables", "PartialProblem", function(object) { variables(object@args[[1]]) })

#' @describeIn PartialProblem Returns the parameters in the problem.
setMethod("parameters", "PartialProblem", function(object) { parameters(object@args[[1]]) })

#' @describeIn PartialProblem Returns the constants in the problem.
setMethod("constants", "PartialProblem", function(object) { constants(object@args[[1]]) })

#' @describeIn PartialProblem Gives the (sub/super)gradient of the expression wrt each variable. Matrix expressions are vectorized, so the gradient is a matrix. NA indicates variable values unknown or outside domain.
setMethod("grad", "PartialProblem", function(object) {
  # Subgrad of g(y) = min f_0(x,y)
  #                   s.t. f_i(x,y) <= 0, i = 1,..,p
  #                        h_i(x,y) == 0, i = 1,...,q
  # Given by Df_0(x^*,y) + \sum_i Df_i(x^*,y) \lambda^*_i
  #          + \sum_i Dh_i(x^*,y) \nu^*_i
  # where x^*, \lambda^*_i, \nu^*_i are optimal primal/dual variables.
  # Add PSD constraints in same way.

  # Short circuit for constant
  if(is_constant(object))
    return(constant_grad(object))

  old_vals <- list()
  for(var in variables(object))
    old_vals[[as.character(id(var))]] <- value(var)

  fix_vars <- list()
  for(var in object@dont_opt_vars) {
    if(is.na(value(var)))
      return(error_grad(object))
    else
      fix_vars <- c(fix_vars, list(var == value(var)))
  }

  prob <- Problem(object@args[[1]]@objective, c(fix_vars, object@args[[1]]@constraints))
  result <- do.call("solve", c(list(a = prob, solver = object@solver), object@.solve_kwargs))

  # Compute gradient.
  if(result$status %in% SOLUTION_PRESENT) {
    sign <- as.numeric(is_convex(object) - is_concave(object))

    # Form Lagrangian.
    lagr <- object@args[[1]]@objective@args[[1]]
    for(constr in object@args[[1]]@constraints) {
      # TODO: better way to get constraint expressions.
      lagr_multiplier <- as.Constant(sign * result$getDualValue(constr))
      lprod <- t(lagr_multiplier) %*% constr@expr
      if(is_scalar(lprod))
        lagr <- lagr + sum(lprod)
      else
        lagr <- lagr + matrix_trace(lprod)
    }

    grad_map <- grad(lagr)   # TODO: After finishing grad implementation, we probably need to update this call by passing in result of solving problem.
    res_grad <- list()
    for(var in object@dont_opt_vars) {
      var_id_char <- as.character(id(var))
      res_grad[[var_id_char]] <- grad_map[[var_id_char]]
    }
  } else   # Unbounded, infeasible, or solver error.
    res_grad <- error_grad(object)

  # Restore the original values to the variables.
  for(var in variables(object))
    value(var) <- old_vals[[as.character(id(var))]]
  return(res_grad)
})

#' @describeIn PartialProblem A list of constraints describing the closure of the region where the expression is finite.
setMethod("domain", "PartialProblem", function(object) {
  # Variables optimized over are replaced in object@args[[1]]
  obj_expr <- object@args[[1]]@objective@args[[1]]
  c(object@args[[1]]@constraints, domain(obj_expr))
})

#' @describeIn PartialProblem Returns the numeric value of the expression.
setMethod("value", "PartialProblem", function(object) {
  old_vals <- list()
  for(var in variables(object))
    old_vals[[as.character(id(var))]] <- value(var)

  fix_vars <- list()
  for(var in object@dont_opt_vars) {
    if(is.na(value(var)))
      return(NA_real_)
    else
      fix_vars <- c(fix_vars, list(var == value(var)))
  }

  prob <- Problem(object@args[[1]]@objective, c(fix_vars, object@args[[1]]@constraints))
  result <- do.call("solve", c(list(a = prob, solver = object@solver), object@.solve_kwargs))

  # Restore the original values to the variables.
  for(var in variables(object))
    value(var) <- old_vals[[as.character(id(var))]]

  # Need to get value returned by solver in case of stacked partial_optimizes.
  return(result$value)
})

#' @describeIn PartialProblem Returns the graph implementation of the object. Chain ids of all the opt_vars.
setMethod("canonicalize", "PartialProblem", function(object) {
  # Canonical form for objective and problem switches from minimize to maximize.
  canon <- canonical_form(object@args[[1]]@objective@args[[1]])
  obj <- canon[[1]]
  constrs <- canon[[2]]
  for(cons in object@args[[1]]@constraints)
    constrs <- c(constrs, list(canonical_form(cons)[[2]]))
  list(obj, constrs)
})


#'
#' The XpressProblem class.
#'
#' This class represents a convex optimization problem associated with the Xpress Optimizer.
#'
#' @slot objective A \linkS4class{Minimize} or \linkS4class{Maximize} object representing the optimization objective.
#' @slot constraints A list of \linkS4class{Constraint} objects representing constraints on the optimization variables.
#' @name XpressProblem-class
#' @aliases XpressProblem
#' @rdname XpressProblem-class
.XpressProblem <- setClass("XpressProblem", representation(.iis = "numeric", .transferRow = "numeric"),
                           prototype(.iis = NA_real_, .transferRow = NA_real_), contains = "Problem")

#' @param objective A \linkS4class{Minimize} or \linkS4class{Maximize} object representing the optimization objective.
#' @param constraints (Optional) A list of \linkS4class{Constraint} objects representing constraints on the optimization variables.
#' @rdname XpressProblem-class
XpressProblem <- function(objective, constraints = list()) {
  .XpressProblem(objective = objective, constraints = constraints)
}

setMethod(".reset_iis", "XpressProblem", function(object) {
  object@.iis <- NA_real_
  object@.transferRow <- NA_real_
  return(object)
})

setMethod("show", "XpressProblem", function(object) {
  cat("XpressProblem(", name(object@objective), ", (", paste(sapply(object@constraints, name), collapse = ", "), "))", sep = "")
})

#'
#' Arithmetic Operations on XpressProblems
#'
#' Add, subtract, multiply, or divide DCP optimization problems.
#'
#' @param e1 The left-hand \linkS4class{XpressProblem} object.
#' @param e2 The right-hand \linkS4class{XpressProblem} object.
#' @return A \linkS4class{XpressProblem} object.
#' @name XpressProblem-arith
NULL

#' @rdname Problem-arith
setMethod("+", signature(e1 = "XpressProblem", e2 = "missing"), function(e1, e2) { XpressProblem(objective = e1@objective, constraints = e1@constraints) })

#' @rdname XpressProblem-arith
setMethod("-", signature(e1 = "XpressProblem", e2 = "missing"), function(e1, e2) { XpressProblem(objective = -e1@objective, constraints = e1@constraints) })

#' @rdname XpressProblem-arith
setMethod("*", signature(e1 = "XpressProblem", e2 = "numeric"), function(e1, e2) {
  XpressProblem(objective = e1@objective * e2, constraints = e1@constraints)
})

#' @rdname XpressProblem-arith
setMethod("+", signature(e1 = "XpressProblem", e2 = "numeric"), function(e1, e2) { if(length(e2) == 1 && e2 == 0) e1 else stop("Unimplemented") })

#' @rdname XpressProblem-arith
setMethod("+", signature(e1 = "numeric", e2 = "XpressProblem"), function(e1, e2) { e2 + e1 })

#' @rdname XpressProblem-arith
setMethod("+", signature(e1 = "XpressProblem", e2 = "XpressProblem"), function(e1, e2) {
  XpressProblem(objective = e1@objective + e2@objective, constraints = unique(c(e1@constraints, e2@constraints)))
})

#' @rdname XpressProblem-arith
setMethod("-", signature(e1 = "XpressProblem", e2 = "numeric"), function(e1, e2) { e1 + (-e2) })

#' @rdname XpressProblem-arith
setMethod("-", signature(e1 = "numeric", e2 = "XpressProblem"), function(e1, e2) { if(length(e1) == 1 && e1 == 0) -e2 else stop("Unimplemented") })

#' @rdname XpressProblem-arith
setMethod("-", signature(e1 = "XpressProblem", e2 = "XpressProblem"), function(e1, e2) {
  XpressProblem(objective = e1@objective - e2@objective, constraints = unique(c(e1@constraints, e2@constraints)))
})

#' @rdname XpressProblem-arith
setMethod("*", signature(e1 = "numeric", e2 = "XpressProblem"), function(e1, e2) { e2 * e1 })

#' @rdname XpressProblem-arith
setMethod("/", signature(e1 = "Problem", e2 = "numeric"), function(e1, e2) {
  XpressProblem(objective = e1@objective * (1.0/e2), constraints = e1@constraints)
})
