#'
#' The Minimize class.
#'
#' This class represents a minimization problem.
#'
#' @slot expr The expression to minimize.
#' @aliases Minimize
#' @export
Minimize <- setClass("Minimize", representation(expr = "ConstValORExpr"), contains = "Canonical")

setMethod("initialize", "Minimize", function(.Object, expr) {
    .Object@expr <- as.Constant(expr)
    if(!all(size(.Object@expr) == c(1,1)))
      stop("The objective must resolve to a scalar")
    return(.Object)
})

setMethod("get_data", "Minimize", function(object) { list() })
setMethod("canonical_form", "Minimize", function(object) { canonicalize(object) })
setMethod("canonicalize", "Minimize", function(object) { canonical_form(object@expr) })
setMethod("variables", "Minimize", function(object) { variables(object@expr) })
setMethod("parameters", "Minimize", function(object) { parameters(object@expr) })
setMethod("constants", "Minimize", function(object) { constants(object@expr) })
setMethod("is_dcp", "Minimize", function(object) { is_convex(object@expr) })
setMethod("value", "Minimize", function(object) { value(object@expr) })
Minimize.primal_to_result <- function(result) { result }

#'
#' The Maximize class.
#'
#' This class represents a maximization problem.
#'
#' @slot expr The expression to maximize.
#' @aliases Maximize
#' @export
Maximize <- setClass("Maximize", contains = "Minimize")

setMethod("-", signature(e1 = "Minimize", e2 = "missing"), function(e1, e2) { Maximize(expr = -e1@expr) })
setMethod("+", signature(e1 = "Minimize", e2 = "Minimize"), function(e1, e2) { Minimize(e1@expr + e2@expr) })
setMethod("+", signature(e1 = "Minimize", e2 = "Maximize"), function(e1, e2) { stop("Problem does not follow DCP rules") })
setMethod("+", signature(e1 = "Minimize", e2 = "numeric"), function(e1, e2) { if(length(e2) == 1 && e2 == 0) e1 else stop("Unimplemented") })
setMethod("+", signature(e1 = "numeric", e2 = "Minimize"), function(e1, e2) { e2 + e1 })
setMethod("-", signature(e1 = "Minimize", e2 = "Minimize"), function(e1, e2) { e1 + (-e2) })
setMethod("-", signature(e1 = "Minimize", e2 = "Maximize"), function(e1, e2) { e1 + (-e2) })
setMethod("-", signature(e1 = "Minimize", e2 = "numeric"), function(e1, e2) { e1 + (-e2) })
setMethod("-", signature(e1 = "numeric", e2 = "Minimize"), function(e1, e2) { if(length(e2) == 1 && e2 == 0) -e1 else stop("Unimplemented") })
setMethod("*", signature(e1 = "Minimize", e2 = "numeric"), function(e1, e2) {
  if(e2 >= 0) Minimize(expr = e1@expr * e2) else Maximize(expr = e1@expr * e2)
})
setMethod("*", signature(e1 = "Maximize", e2 = "numeric"), function(e1, e2) {
  if(e2 < 0) Minimize(expr = e1@expr * e2) else Maximize(expr = e1@expr * e2)
})
setMethod("*", signature(e1 = "numeric", e2 = "Minimize"), function(e1, e2) { e2 * e1 })
setMethod("/", signature(e1 = "Minimize", e2 = "numeric"), function(e1, e2) { e1 * (1.0/e2) })

setMethod("-", signature(e1 = "Maximize", e2 = "missing"), function(e1, e2) { Minimize(expr = -e1@expr) })
setMethod("+", signature(e1 = "Maximize", e2 = "Maximize"), function(e1, e2) { Maximize(expr = e1@expr + e2@expr) })
setMethod("+", signature(e1 = "Maximize", e2 = "Minimize"), function(e1, e2) { stop("Problem does not follow DCP rules") })

setMethod("canonicalize", "Maximize", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constraints <- canon[[2]]
  list(neg_expr(obj), constraints)
})

setMethod("is_dcp", "Maximize", function(object) { is_concave(object@expr) })
setMethod("is_quadratic", "Maximize", function(object) { is_quadratic(object@expr) })
Maximize.primal_to_result <- function(result) { -result }

.SolverStats <- setClass("SolverStats", representation(solver_name = "character", solve_time = "numeric", setup_time = "numeric", num_iters = "numeric"), 
                         prototype(solver_name = NA_character_, solve_time = NA_real_, setup_time = NA_real_, num_iters = NA_real_))
SolverStats <- function(results_dict = list(), solver_name = NA_character_) {
  solve_time <- NA_real_
  setup_time <- NA_real_
  num_iters <- NA_real_
  
  if(SOLVE_TIME %in% names(results_dict))
    solve_time <- results_dict[SOLVE_TIME]
  if(SETUP_TIME %in% names(results_dict))
    setup_time <- results_dict[SETUP_TIME]
  if(NUM_ITERS %in% names(results_dict))
    num_iters <- results_dict[NUM_ITERS]
  .SolverStats(solver_name = solver_name, solve_time = solve_time, setup_time = setup_time, num_iters = num_iters)
}

.SizeMetrics <- setClass("SizeMetrics", representation(num_scalar_variables = "numeric", num_scalar_data = "numeric", num_scalar_eq_constr = "numeric", num_scalar_leq_constr = "numeric",
                                                       max_data_dimension = "numeric", max_big_small_squared = "numeric"),
                         prototype(num_scalar_variables = NA_real_, num_scalar_data = NA_real_, num_scalar_eq_constr = NA_real_, num_scalar_leq_constr = NA_real_,
                                   max_data_dimension = NA_real_, max_big_small_squared = NA_real_))

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
    if(is(constraint, "LeqConstraint"))
      num_scalar_leq_constr <- num_scalar_leq_constr + prod(size(constraint@.expr))
  }
  
  .SizeMetrics(num_scalar_variables = num_scalar_variables, num_scalar_data = num_scalar_data, num_scalar_eq_constr = num_scalar_eq_constr, num_scalar_leq_constr = num_scalar_leq_constr,
               max_data_dimension = max_data_dimension, max_big_small_squared = max_big_small_squared)
}

setClassUnion("SizeMetricsORNull", c("SizeMetrics", "NULL"))

#'
#' The Problem class.
#'
#' This class represents the convex optimization problem.
#'
#' @slot objective The expression to minimize or maximize.
#' @slot constraints (Optional) A list of constraints on the problem variables.
#' @aliases Problem
#' @export
.Problem <- setClass("Problem", representation(objective = "Minimize", constraints = "list", value = "numeric", status = "character", .cached_data = "list", .separable_problems = "list", .size_metrics = "SizeMetricsORNull", .solver_stats = "list"),
                    prototype(constraints = list(), value = NA_real_, status = NA_character_, .cached_data = list(), .separable_problems = list(), .size_metrics = NULL, .solver_stats = list()),
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

Problem <- function(objective, constraints = list()) {
  .Problem(objective = objective, constraints = constraints)
}

CachedProblem <- function(objective, constraints) { list(objective = objective, constraints = constraints, class = "CachedProblem") }
SolveResult <- function(opt_value, status, primal_values, dual_values) { list(opt_value = opt_value, status = status, primal_values = primal_values, dual_values = dual_values, class = "SolveResult") }

setMethod("initialize", "Problem", function(.Object, ..., objective, constraints = list(), value = NA_real_, status = NA_character_, .cached_data = list(), .separable_problems = list(), .size_metrics = SizeMetrics(), .solver_stats = list()) {
  .Object@objective <- objective
  .Object@constraints <- constraints
  .Object@value <- value
  .Object@status <- status

  # Cached processed data for each solver.
  .reset_cache <- function(object) {
    for(solver_name in sapply(SOLVERS, name))
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

setMethod("value", "Problem", function(object) { object@value })
setMethod("status", "Problem", function(object) { object@status })
setMethod("is_dcp", "Problem", function(object) {
  all(sapply(c(object@constraints, list(object@objective)), is_dcp))
})

setMethod("canonicalize", "Problem", function(object) {
  obj_canon <- canonical_form(object@objective)
  canon_constr <- obj_canon[[2]]
  
  for(constr in object@constraints)
    canon_constr <- c(canon_constr, canonical_form(constr)[[2]])
  list(obj_canon[[1]], canon_constr)
})

setMethod("variables", "Problem", function(object) {
  vars_ <- variables(object@objective)
  constrs_ <- lapply(object@constraints, function(constr) { variables(constr) })
  unique(flatten_list(c(vars_, constrs_)))
})

setMethod("parameters", "Problem", function(object) {
  params <- parameters(object@objective)
  constrs_ <- lapply(object@constraints, function(constr) { parameters(constr) })
  unique(flatten_list(c(params, constrs_)))
})

setMethod("constants", "Problem", function(object) {
  constants_ <- lapply(object@constraints, function(constr) { constants(constr) })
  constants_ <- c(constants(object@objective), constants_)
  unique(flatten_list(constants_))   # TODO: Check duplicated constants are removed correctly
})

setMethod("size_metrics", "Problem", function(object) { object@.size_metrics })
setMethod("solver_stats", "Problem", function(object) { object@.solver_stats })

# solve.Problem <- function(object, method, ...) {
#  if(missing(method))
#    .solve(object, ...)
#  else {
#    func <- Problem.REGISTERED_SOLVE_METHODS[func_name]
#    func(object, ...)
#  }
# }

# Problem.register_solve <- function(cls, name, func) {
#  cls.REGISTERED_SOLVE_METHODS[name] <- func
# }

setMethod("get_problem_data", signature(object = "Problem", solver = "character"), function(object, solver) {
  canon <- canonicalize(object)
  objective <- canon[[1]]
  constraints <- canon[[2]]
  
  # Raise an error if the solver cannot handle the problem
  validate_solver(SOLVERS[[solver]], constraints)
  Solver.get_problem_data(SOLVERS[[solver]], objective, constraints, object@.cached_data)
})

solve.Problem <- function(object, solver = ECOS_NAME, ignore_dcp = FALSE, warm_start = FALSE, verbose = FALSE, parallel = FALSE, ...) {
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
    if(objective != object@.cached_data[PARALLEL]@objective ||
       constraints != object@.cached_data[PARALLEL]@constraints) {
      object@.separable_problems <- get_separable_problems(object)
      object@.cached_data[PARALLEL] <- CachedProblem(objective, constraints)
    }
    
    if(length(object@.separable_problems) > 1)
      .parallel_solve(object, solver, ignore_dcp, warm_start, verbose, ...)
  }
  # sym_data <- get_sym_data(solver, objective, constraints, object@.cached_data)
  
  print("Calling CVXcanon")
  if(is(object@objective, "Minimize")) {
    sense <- "Minimize"
    canon_objective <- objective
  } else {
    sense <- "Maximize"
    canon_objective <- neg_expr(objective)  # preserve sense
  }
  
  solve_int(sense, canon_objective, constraints, verbose, ...)
}

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
          subproblem <- .save_value(var, primal_value)   # TODO: Fix this since R makes copies
        }
        
        for(j in 1:length(solve_result$dual_values)) {
          constr <- subproblem@constraints[j]
          dual_value <- solve_result$dual_values[j]
          subproblem <- .save_value(constr, dual_value)   # TODO: Fix this since R makes copies
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

.update_problem_state.Problem <- function(object, results_dict, sym_data, solver) {
  if(results_dict[STATUS] %in% SOLUTION_PRESENT) {
    object <- .save_values(results_dict[PRIMAL], variables(object), sym_data@var_offsets)
    # Not all solvers provide dual variables
    if(EQ_DUAL %in% names(results_dict))
      object <- .save_dual_values(object, results_dict[EQ_DUAL], sym_data@constr_map[EQ], list(EqConstraint))
    if(INEQ_DUAL %in% names(results_dict))
      object <- .save_dual_values(object, results_dict[INEQ_DUAL], sym_data@constr_map[LEQ], list(LeqConstraints, PSDConstraint))
    
    # Correct optimal value if the objective was Maximize
    value <- results_dict[VALUE]
    object@value <- primal_to_result(object@objective, value)
  } else if(results_dict[STATUS] %in% INF_OR_UNB)   # Infeasible or unbounded
    object <- .handle_no_solution(object, results_dict[STATUS])
  else   # Solver failed to solve
    stop("Solver failed. Try another.")
  object@status <- results_dict[STATUS]
  object@solver_stats <- SolverStats(results_dict, solver)
  object
}

unpack_results.Problem <- function(object, solve, results_dict) {
  canon <- canonicalize(object)
  objective <- canon[[1]]
  constraints <- canon[[2]]
  sym_data <- get_sym_data(solver, objective, constraints, object@.cached_data)
  data <- list(sym_data@dims, 0)
  names(data) <- c(DIMS, OFFSET)
  results_dict <- format_results(solver, results_dict, data, object@.cached_data)
  .update_problem_state(object, results_dict, sym_data, solver)
}

.handle_no_solution.Problem <- function(object, status) {
  # Set all primal and dual variable values to NA
  # TODO: This won't work since R creates copies
  for(var_ in variables(object))
    object <- .save_value(var_, NA)
  
  for(i in 1:length(object@constraints))
    object@constraints[i] <- .save_value(object@constraints[i], NA)
  
  # Set the problem value
  if(status %in% c("INFEASIBLE", "INFEASIBLE_INACCURATE"))
    object@value <- primal_to_result(object@objective, Inf)
  else if(status %in% c("UNBOUNDED", "UNBOUNDED_INACCURATE"))
    object@value <- primal_to_result(object@objective, -Inf)
  object
}

.save_dual_values.Problem <- function(object, result_vec, constraints, constr_types) {
  constr_offsets <- list()
  offset <- 0
  for(constr in constraints) {
    constr_offsets[constr@constr_id] <- offset
    offset <- offset + prod(size(constr))
  }
  
  active_constraints <- list()
  for(constr in object@constraints) {
    # Ignore constraints of the wrong type
    if(class(constr) %in% constr_types)
      active_constraints <- c(active_constraints, constr)
  }
  .save_values(object, result_vec, active_constraints, constr_offsets)
}

.save_values.Problem <- function(object, result_vec, objects, offset_map) {
  if(length(result_vec) > 0)   # Cast to desired matrix type
    result_vec <- as.matrix(result_vec)
  
  objects_new <- list()
  for(obj in objects) {
    size <- size(obj)
    rows <- size[1]
    cols <- size[2]
    id_char <- as.character(obj@id)
    if(id_char %in% names(offset_map)) {
      offset <- offset_map[id_char]
      
      # Handle scalars
      if(all(c(rows, cols) == c(1,1)))
        value <- index(result_vec, c(offset, 0))
      else {
        value <- matrix(0, nrow = rows, ncol = cols)
        value <- block_add(value, result_vec[offset:(offset + rows*cols)], 1, 1, rows, cols)
        offset <- offset + rows * cols
      }
    } else  # The variable was multiplied by zero
      value <- matrix(0, nrow = rows, ncol = cols)
    obj <- save_value(obj, value)
    objects_new <- c(objects_new, obj)
  }
  objects_new
}

setMethod("+", signature(e1 = "Problem", e2 = "missing"), function(e1, e2) { Problem(objective = e1@objective, constraints = e1@constraints) })
setMethod("-", signature(e1 = "Problem", e2 = "missing"), function(e1, e2) { Problem(objective = -e1@objective, constraints = e1@constraints) })
setMethod("+", signature(e1 = "Problem", e2 = "numeric"), function(e1, e2) { if(length(e2) == 1 && e2 == 0) e1 else stop("Unimplemented") })
setMethod("+", signature(e1 = "numeric", e2 = "Problem"), function(e1, e2) { e2 + e1 })
setMethod("+", signature(e1 = "Problem", e2 = "Problem"), function(e1, e2) {
  Problem(objective = e1@objective + e2@objective, constraints = unique(c(e1@constraints, e2@constraints)))
})
setMethod("-", signature(e1 = "Problem", e2 = "numeric"), function(e1, e2) { e1 + (-e2) })
setMethod("-", signature(e1 = "numeric", e2 = "Problem"), function(e1, e2) { if(length(e1) == 1 && e2 == 0) -e2 else stop("Unimplemented") })
setMethod("-", signature(e1 = "Problem", e2 = "Problem"), function(e1, e2) {
  Problem(objective = e1@objective - e2@objective, constraints = unique(c(e1@constraints, e2@constraints)))
})
setMethod("*", signature(e1 = "Problem", e2 = "numeric"), function(e1, e2) {
  Problem(objective = e1@objective * e2, constraints = e1@constraints)
})
setMethod("*", signature(e1 = "numeric", e2 = "Problem"), function(e1, e2) { e2 * e1 })
setMethod("/", signature(e1 = "Problem", e2 = "numeric"), function(e1, e2) {
  Problem(objective = e1@objective * (1.0/e2), constraints = e1@constraints)
})
