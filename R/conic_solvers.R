is_stuffed_cone_constraint <- function(constraint) {
  # Conic solvers require constraints to be stuffed in the following way.
  if(length(variables(constraint)) != 1)
    return(FALSE)
  for(arg in constraints@args) {
    if(class(arg) == "Reshape")
      arg <- arg@args[[1]]
    if(class(arg) == "AddExpression") {
      if(!class(arg@args[[1]]) %in% c("MulExpression", "Multiply"))
        return(FALSE)
      if(class(arg@args[[1]]@args[[1]]) != "Constant")
        return(FALSE)
      if(class(arg@args[[2]]) != "Constant")
        return(FALSE)
    } else if(class(arg) %in% c("MulExpression", "Multiply")) {
      if(class(arg@args[[1]]) != "Constant")
        return(FALSE)
    } else
      return(FALSE)
  }
  return(TRUE)
}

is_stuffed_cone_objective <- function(objective) {
  # Conic solvers require objectives to be stuffed in the following way.
  expr <- objective@expr
  return(is_affine(expr) && length(variables(expr)) == 1 && class(expr) == "AddExpression" && length(expr@args) == 2
                         && class(expr@args[[1]]) %in% c("MulExpression", "Multiply") && class(expr@args[[2]]) == "Constant")
}


#' Summary of cone dimensions present in constraints.
#'
#'    Constraints must be formatted as dictionary that maps from
#'    constraint type to a list of constraints of that type.
#'
#'    Attributes
#'    ----------
#'    zero : int
#'        The dimension of the zero cone.
#'    nonpos : int
#'        The dimension of the non-positive cone.
#'    exp : int
#'        The dimension of the exponential cone.
#'    soc : list of int
#'        A list of the second-order cone dimensions.
#'    psd : list of int
#'        A list of the positive semidefinite cone dimensions, where the
#'        dimension of the PSD cone of k by k matrices is k.
.ConeDims <- setClass("ConeDims", representation(constr_map = "list", zero = "numeric", nonpos = "numeric", exp = "numeric", soc = "numeric", psd = "numeric"),
                                  prototype(zero = NA_real_, nonpos = NA_real_, exp = NA_real_, soc = NA_real_, psd = NA_real_))

ConeDims <- function(constr_map) { .ConeDims(constr_map = constr_map) }

setMethod("initialize", "ConeDims", function(.Object, constr_map, zero = NA_real_, nonpos = NA_real_, exp = NA_real_, soc = NA_real_, psd = NA_real_) {
  .Object@zero <- sum(sapply(constr_map$Zero, function(c) { size(c) }))
  .Object@nonpos <- sum(sapply(constr_map$NonPos, function(c) { size(c) }))
  .Object@exp <- sum(sapply(constr_map$ExpCone, function(c) { num_cones(c) }))
  .Object@soc <- sapply(constr_map$SOC, function(c) { cone_sizes(c) })
  .Object@psd <- sapply(constr_map$PSD, function(c) { dim(c)[1] })
  return(.Object)
})

# Conic solver class with reduction semantics.
ConicSolver <- setClass("ConicSolver", contains = "Solver")

# The key that maps to ConeDims in the data returned by apply().
setMethod("dims", "ConicSolver", function(object) { "dims" })

# Every conic solver must support Zero and NonPos constraints.
setMethod("supported_constraints", "ConicSolver", function(object) { c("Zero", "NonPos") })

# Some solvers cannot solve problems that do not have constraints.
# For such solvers, requires_constr should return TRUE.
setMethod("requires_constr", "ConicSolver", function(object) { FALSE })

setMethod("accepts", signature(object = "ConicSolver", problem = "Problem"), function(object, problem) {
  return(class(problem@objective) == "Minimize" && (mip_capable(object) || !is_mixed_integer(problem)) && is_stuffed_cone_objective(problem@objective)
    && !convex_attributes(variables(problem)) && (length(problem@constraints) > 0 || !requires_constr(object))
    && all(sapply(problem@constraints, function(c) { class(c) %in% supported_constraints(object) }))
    && all(sapply(problem@constraints, function(c) { is_stuffed_cone_constraints(c) })))
})

ConicSolver.get_coeff_offset <- function(expr) {
  # Return the coefficient and offset in A %*% x + b.
  if(class(expr) == "Reshape")   # May be a Reshape as root.
    expr <- expr@args[[1]]
  if(length(expr@args[[1]]@args) == 0) {   # Convert data to float64.
    # expr is t(c) %*% x
    offset <- 0
    coeff <- as.numeric(value(expr@args[[1]]))
  } else {
    # expr is t(c) %*% x + d
    offset <- as.numeric(value(expr@args[[2]]))
    coeff <- as.numeric(value(expr@args[[1]]@args[[1]]))
  }
  # Convert scalars to sparse matrices.
  if(is.atomic(coeff) && length(coeff) == 1)
    coeff <- Matrix(coeff, sparse = TRUE)
  return(list(coeff, offset))
}

ConicSolver.get_spacing_matrix <- function(shape, spacing, offset) {
  # Returns a sparse matrix that spaces out an expression.
  val_arr <- c()
  row_arr <- c()
  col_arr <- c()

  # Selects from each column.
  for(var_row in 1:shape[2]) {
    val_arr <- c(val_arr, 1.0)
    row_arr <- c(row_arr, spacing*var_row + offset)
    col_arr <- c(col_arr, var_row)
  }
  return(sparseMatrix(i = row_arr, j = col_arr, x = val_arr))
}

setMethod("format_constr", "ConicSolver", function(object, problem, constr, exp_cone_order) {
  coeffs <- list()
  offsets <- list()
  for(arg in constr@args) {
    res <- ConicSolver.get_coeff_offset(arg)
    coeffs <- c(coeffs, Matrix(res[[1]], sparse = TRUE))
    offsets <- c(offsets, res[[2]])
  }
  height <- sum(sapply(coeffs, function(c) { dim(c)[1] }))

  if(class(constr) %in% c("NonPos", "Zero"))
    # Both of these constraints have but a single argument.
    # t(c) %*% x + b (<)= 0 if and only if t(c) %*% x (<)= b.
    return(list(coeffs[[1]], -offsets[1]))
  else if(class(constr) == "SOC") {
    # Group each t row with appropriate X rows.
    mat_arr <- list()
    offset <- rep(0, height)
    if(constr@axis == 2)
      gap <- nrow(constr@args[[2]]) + 1
    else
      gap <- ncol(constr@args[[2]]) + 1

    for(i in 1:size(constr@args[[1]])) {
      offset[(i-1)*gap+1] <- offsets[1][i]
      c(mat_arr, coeffs[1][i,])
      if(constr@axis == 2) {
        offset[((i-1)*gap+2):(i*gap+1)] <- offsets[2][((i-1)*(gap-1)+1):(i*(gap-1)+1),]
        mat_arr <- c(mat_arr, coeffs[2][((i-1)*(gap-1)+1):(i*(gap-1)+1),])
      } else {
        offset[((i-1)*gap+2):(i*gap+1)] <- offsets[2][seq(i, length(offsets[2]), gap-1)]
        mat_arr <- c(mat_arr, coeffs[2][seq(i, nrow(coeffs[2]), gap-1),])
      }
    }
    return(list(-Matrix(do.call(rbind, mat_arr), sparse = TRUE), offset))
  } else if(class(constr) == "ExpCone") {
    for(i in 1:length(coeffs)) {
      mat <- ConicSolver.get_spacing_matrix(c(height, nrow(coeffs[[i]])), length(exp_cone_order), exp_cone_order[i])
      offsets[i] <- mat %*% offsets[i]
      coeffs[i] <- -mat %*% coeffs[i]
    }
    # return(list(sum(coeffs), sum(offsets)))
    return(list(Reduce("+", coeffs), Reduce("+", offsets)))
  } else
    # subclasses must handle PSD constraints.
    stop("Unsupported constraint type.")
})

setMethod("group_coeff_offset", "ConicSolver", function(object, problem, constraints, exp_cone_order) {
  # Combine the constraints into a single matrix, offset.
  if(is.na(constraints) || is.null(constraints) || length(constraints) == 0)
    return(list(NA, NA))

  matrices <- list()
  offsets <- list()
  for(cons in constraints) {
    res <- format_constr(object, problem, cons, exp_cone_order)
    matrices <- c(matrices, res[[1]])
    offsets <- c(offsets, res[[2]])
  }
  coeff <- Matrix(do.call(rbind, matrices), sparse = TRUE)
  offset <- do.call(cbind, offsets)
  return(list(coeff, offset))
})

setMethod("invert", "ConicSolver", signature(object = "ConicSolver", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  # Returns the solution to the original problem given the inverse_data.
  status <- solution$status
  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution$value
    primal_vars <- list()
    primal_vars[inverse_data[var_id(object)]] <- solution$primal
    eq_dual <- get_dual_values(solution$eq_dual, extract_dual_value, inverse_data[eq_constr(object)])
    leq_dual <- get_dual_values(solution$ineq_dual, extract_dual_value, inverse_data[neq_constr(object)])
    eq_dual <- modifyList(eq_dual, leq_dual)
    dual_vars <- eq_dual
  } else {
    if(status == INFEASIBLE)
      opt_val <- Inf
    else if(status == UNBOUNDED)
      opt_val <- -Inf
    else
      opt_val <- NA
    primal_vars <- NA
    dual_vars <- NA
  }
  return(Solution(status, opt_val, primal_vars, dual_vars, list()))
})

CBC <- setClass("CBC", contains = "ConicSolver")

# Solver capabilities.
setMethod("mip_capable", "CBC", function(object) { TRUE })

# Map of GLPK MIP status to CVXR status.
setMethod("status_map_mip", "CBC", function(object) {
  list(solution = OPTIMAL, relaxation_infeasible = INFEASIBLE, stopped_on_user_event = SOLVER_ERROR)
})

setMethod("status_map_lp", "CBC", function(object) {
  list(optimal = OPTIMAL, primal_infeasible = INFEASIBLE, stopped_due_to_errors = SOLVER_ERROR,
       stopped_by_event_handler = SOLVER_ERROR)
})

setMethod("name", "CBC", function(x) { CBC_NAME })
setMethod("import_solver", "CBC", function(object) {
  requireNamespace("rcbc", quietly = TRUE)
})

setMethod("accepts", signature(object = "CBC", problem = "Problem"), function(object, problem) {
  # Can CBC solve the problem?
  # TODO: Check if the matrix is stuffed.
  if(!is_affine(problem@objective@args[[1]]))
    return(FALSE)
  for(constr in problem@constraints) {
    if(!class(constr) %in% supported_constraints(object))
      return(FALSE)
    for(arg in constr@args) {
      if(!is_affine(arg))
        return(FALSE)
    }
  }
  return(TRUE)
})

setMethod("apply", signature(object = "CBC", problem = "Problem"), function(object, problem) {
  # Returns a new problem and data for inverting the new solution.
  data <- list()
  objective <- canonical_form(problem@objective)[[1]]
  constraints <- lapply(problem@constraints, function(c) { canonical_form[[2]] })
  constraints <- unlist(constraints, recursive = TRUE)
  data$objective <- objective
  data$constraints <- constraints
  variables <- variables(problem)[[1]]
  data[BOOL_IDX] <- lapply(variables@boolean_idx, function(t) { as.integer(t[[1]]) })
  data[INT_IDX] <- lapply(variables@integer_idx, function(t) { as.integer(t[[1]]) })

  # Order and group constraints.
  inv_data <- list()
  inv_data[var_id(object)] <- id(variables(problem)[[1]])
  eq_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "Zero" })]
  inv_data[eq_constr(object)] <- eq_constr
  leq_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "NonPos" })]
  inv_data[neq_constr(object)] <- leq_constr
  return(list(data, inv_data))
})

setMethod("invert", signature(object = "CBC", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  # Returns the solution to the original problem given the inverse_data.
  status <- solution$status

  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution$value
    primal_vars[var_id(object)] <- solution$primal
  } else {
    if(status == INFEASIBLE)
      opt_val <- Inf
    else if(status == UNBOUNDED)
      opt_val <- -Inf
    else
      opt_val <- NA
    primal_vars <- NA
  }
  dual_vars <- NA

  return(Solution(status, opt_val, primal_vars, dual_vars, list()))
})

setMethod("solve_via_data", "CBC", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  solver <- CBC_OLD()
  solver_opts[BOOL_IDX] <- data[BOOL_IDX]
  solver_opts[INT_IDX] <- data[INT_IDX]
  prob_data <- list()
  prob_data[name(object)] <- ProblemData()
  return(solve(solver, data$objective, data$constraints, prob_data, warm_start, verbose, solver_opts))
})

CPLEX <- setClass("CPLEX", contains = "ConicSolver")

setMethod("mip_capable", "CPLEX", function(object) { TRUE })
setMethod("supported_constraints", "CPLEX", function(object) { c(supported_constraints(ConicSolver()), "SOC") })
setMethod("name", "CPLEX", function(x) { CPLEX_NAME })
setMethod("import_solver", "CPLEX", function(object) { requireNamespace("Rcplex", quietly = TRUE) })
setMethod("accepts", signature(object = "CPLEX", problem = "Problem"), function(object, problem) {
  # Can CPLEX solve the problem?
  # TODO: Check if the matrix is stuffed.
  if(!is_affine(problem@objective@args[[1]]))
    return(FALSE)
  for(constr in problem@constraints) {
    if(!class(constr) %in% supported_constraints(object))
      return(FALSE)
    for(arg in constr@args) {
      if(!is_affine(arg))
        return(FALSE)
    }
  }
  return(TRUE)
})

setMethod("apply", signature(object = "CPLEX", problem = "Problem"), function(object, problem) {
  # Returns a new problem and data for inverting the new solution.
  data <- list()
  objective <- canonical_form(problem@objective)[[1]]
  constraints <- lapply(problem@constraints, function(c) { canonical_form(c)[[2]] })
  constraints <- unlist(constraints, recursive = TRUE)
  data$objective <- objective
  data$constraints <- constraints
  variables <- variables(problem)[[1]]
  data[BOOL_IDX] <- lapply(variables@boolean_idx, function(t) { t[[1]] })
  data[INT_IDX] <- lapply(variables@integer_idx, function(t) { t[[1]] })

  # Order and group constraints.
  inv_data <- list()
  inv_data[var_id(object)] <- id(variables(problem)[[1]])
  eq_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "Zero" })]
  inv_data[eq_constr(object)] <- eq_constr
  leq_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "NonPos" })]
  soc_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "SOC" })]
  inv_data[neq_constr(object)] <- c(leq_constr, soc_constr)
  inv_data$is_mip <- length(data[BOOL_IDX]) > 0 || length(data[INT_IDX]) > 0
  return(list(data, inv_data))
})

setMethod("invert", signature(object = "CPLEX", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  # Returns the solution to the original problem given the inverse_data.
  status <- solution$status

  primal_vars <- NA
  dual_vars <- NA
  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution$value
    primal_vars <- list()
    primal_vars[inverse_data[var_id(object)]] <- solution$primal
    if(!inverse_data$is_mip) {
      eq_dual <- get_dual_values(solution$eq_dual, extract_dual_value, inverse_data[eq_constr(object)])
      leq_dual <- get_dual_values(solution$ineq_dual, extract_dual_value, inverse_data[neq_constr(object)])
      eq_dual <- modifyList(eq_dual, leq_dual)
      dual_vars <- eq_dual
    } else {
      if(status == INFEASIBLE)
        opt_val <- Inf
      else if(status == UNBOUNDED)
        opt_val <- -Inf
      else
        opt_val <- NA_real_
    }
  }

  return(Solution(status, opt_val, primal_vars, dual_vars, list()))
})

setMethod("solve_via_data", "CPLEX", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  solver <- CPLEX_OLD()
  solver_opts[BOOL_IDX] <- data[BOOL_IDX]
  solver_opts[INT_IDX] <- data[INT_IDX]
  prob_data <- list()
  prob_data[name(object)] <- ProblemData()
  solve(solver, data$objective, data$constraints, prob_data, warm_start, verbose, solver_opts)
})

CVXOPT <- setClass("CVXOPT", contains = "ConicSolver")

# Solver capabilities.
setMethod("mip_capable", "CVXOPT", function(object) { FALSE })
setMethod("supported_constraints", "CVXOPT", function(object) { c(supported_constraints(ConicSolver()), "SOC", "ExpCone", "PSD") })

# Map of CVXOPT status to CVXR status.
setMethod("status_map", "CVXOPT", function(object, status) {
  list(optimal = OPTIMAL, infeasible = INFEASIBLE, unbounded = UNBOUNDED, solver_error = SOLVER_ERROR)
})

setMethod("name", "CVXOPT", function(x) { CVXOPT_NAME })
setMethod("import_solver", "CVXOPT", function(object) { requireNamespace("cccopt", quietly = TRUE) })

setMethod("accepts", signature(object = "CVXOPT", problem = "Problem"), function(object, problem) {
  # Can CVXOPT solver the problem?
  # TODO: Check if the matrix is stuffed.
  if(!is_affine(problem@objective@args[[1]]))
    return(FALSE)
  for(constr in problem@constraints) {
    if(!class(constr) %in% supported_constraints(object))
      return(FALSE)
    for(arg in constr@args) {
      if(!is_affine(arg))
        return(FALSE)
    }
  }
  return(TRUE)
})

setMethod("apply", signature(object = "CVXOPT", problem = "Problem"), function(object, problem) {
  data <- list()
  objective <- canonical_form(problem@objective)[[1]]
  constraints <- lapply(problem@constraints, function(c) { canonical_form(c)[[1]] })
  constraints <- unlist(constraints, recursive = TRUE)
  data$objective <- objective
  data$constraints <- constraints
  data[dims(ConicSolver())] <- ConeDims(group_constraints(problem@constraints))
  variables <- variables(problem)[[1]]
  data[BOOL_IDX] <- lapply(variables@boolean_idx, function(t) { t[[1]] })
  data[INT_IDX] <- lapply(variables@integer_idx, function(c) { t[[1]] })

  inv_data <- list()
  inv_data[var_id(object)] <- id(variables(problem)[[1]])

  # Order and group constraints.
  eq_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "Zero" })]
  inv_data[eq_constr(CVXOPT())] <- eq_constr
  leq_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "NonPos" })]
  soc_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "SOC" })]
  sdp_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "PSD" })]
  exp_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "ExpCone" })]
  inv_data[neq_constr(CVXOPT())] <- c(leq_constr, soc_constr, sdp_constr, exp_constr)
  return(list(data, inv_data))
})

setMethod("solve_via_data", "CVXOPT", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  solver <- CVXOPT_OLD()
  prob_data <- list()
  prob_data[name(object)] <- ProblemData()
  solve(solver, data$objective, data$constraints, prob_data, warm_start, verbose, solver_opts)
})

ECOS_BB <- setClass("ECOS_BB", contains = "ECOS")

setMethod("mip_capable", "ECOS_BB", function(object) { TRUE })
setMethod("name", "ECOS_BB", function(x) { ECOS_BB_NAME })
setMethod("apply", signature(object = "ECOS_BB", problem = "Problem"), function(object, problem) {
  res <- callNextMethod(object, problem)
  data <- res[[1]]
  inv_data <- res[[2]]

  # Because the problem variable is single dimensional, every
  # boolean/integer index has length one.
  var <- variables(problem)[[1]]
  data[BOOL_IDX] <- lapply(var@boolean_idx, function(t) { as.integer(t[1]) })
  data[INT_IDX] <- lapply(var@integer_idx, function(t) { as.integer(t[1]) })
  return(list(data, inv_data))
})

setMethod("solve_via_data", "ECOS_BB", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  requireNamespace("ECOSolveR", quietly = TRUE)
  cones <- dims_to_solver_dict(data[dims(ConicSolver())])

  # Default verbose to false for BB wrapper.
  if("mi_verbose" %in% names(solver_opts)) {
    mi_verbose <- solver_opts$mi_verbose
    solver_opts$mi_verbose <- NULL
  } else
    mi_verbose <- verbose
  solution <- ECOSolveR::solve(data[C_KEY], data[G_KEY], data[H_KEY], cones, data[A_KEY], data[B_KEY],
                               verbose = verbose, mi_verbose = mi_verbose, bool_vars_idx = data[BOOL_IDX], int_vars_id = data[INT_IDX], solver_opts)
  return(solution)
})

# Utility method for formatting a ConeDims instance into a dictionary
# that can be supplied to ECOS.
dims_to_solver <- function(cone_dims) {
  cones <- list(l = as.integer(cone_dims@nonpos),
                q = lapply(cone_dims@soc, function(v) { as.integer(v) }),
                e = as.integer(cone_dims@exp))
  return(cones)
}

ECOS <- setClass("ECOS", contains = "ConicSolver")

# Solver capabilities.
setMethod("mip_capable", "ECOS", function(object) { FALSE })
setMethod("supported_constraints", "ECOS", function(object) { c(supported_constraints(ConicSolver()), "SOC", "ExpCone") })

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
setMethod("status_map", "ECOS", function(object, status) {
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

# Order of exponential cone arguments for solver.
setMethod("exp_cone_order", "ECOS", function(object) { c(0, 2, 1) })

setMethod("import_solver", "ECOS", function(object) {
  requireNamespace("ECOSolveR", quietly = TRUE)
})

setMethod("name", "ECOS", function(x) { ECOS_NAME })
setMethod("apply", "ECOS", signature(object = "ECOS", problem = "Problem"), function(object, problem) {
  data <- list()
  inv_data <- list()
  inv_data[var_id(object)] <- id(variables(problem)[[1]])
  offsets <- get_coeff_offset(ConicSolver(), problem@objective@args[[1]])
  data[C_KEY] <- as.vector(offsets[[1]])
  data[OFFSET] <- offsets[[2]]
  inv_data[OFFSET] <- data[OFFSET][[1]]

  constr_map <- group_constraints(problem@constraints)
  data[dims(ConicSolver())] <- ConeDims(constr_map)

  inv_data[eq_constr(object)] <- constr_map$Zero
  offsets <- group_coeff_offset(object, problem, constr_map$Zero, exp_cone_order(ECOS()))
  data[A_KEY] <- offsets[[1]]
  data[[B_KEY]] <- offsets[[2]]

  # Order and group nonlinear constraints.
  neq_constr <- c(constr_map$NonPos, constr_map$SOC, constr_map$ExpCone)
  inv_data[neq_constr(object)] <- neq_constr
  offsets <- group_coeff_offset(object, problem, neq_constr, exp_cone_order(ECOS()))
  data[G_KEY] <- offsets[[1]]
  data[H_KEY] <- offsets[[2]]

  return(list(data, inv_data))
})

setMethod("invert", signature(object = "ECOS", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  status <- status_map(object, solution$info$exitFlag)

  # Timing data.
  attr <- list()
  attr[SOLVE_TIME] <- solution$info$timing$tsolve
  attr[SETUP_TIME] <- solution$info$timing$tsetup
  attr[NUM_ITERS] <- solution$info$iter

  if(status %in% SOLUTION_PRESENT) {
    primal_val <- solution$info$pcost
    opt_val <- primal_val + inverse_data[OFFSET]
    primal_vars <- list()
    primal_vars[inverse_data[var_id(object)]] <- as.matrix(solution$x)

    eq_dual <- get_dual_values(solution$y, extract_dual_value, inverse_data[eq_constr(object)])
    leq_dual <- get_dual_values(solution$z, extract_dual_value, inverse_data[neq_constr(object)])
    eq_dual <- modifyList(eq_dual, leq_dual)
    dual_vars <- eq_dual

    return(Solution(status, opt_val, primal_vars, dual_vars, attr))
  } else
    return(failure_solution(status))
})

setMethod("solve_via_data", "ECOS", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  requireNamespace("ECOSolveR", quietly = TRUE)
  cones <- dims_to_solver_dict(data[dims(ConicSolver())])
  solution <- ECOSolveR::solve(data[C_KEY], data[G_KEY], data[H_KEY], cones, data[A_KEY], data[B_KEY], verbose = verbose, solver_opts)
  return(solution)
})

Elemental <- setClass("Elemental", contains = "ConicSolver")

# Solver capabilities
setMethod("mip_capable", "Elemental", function(object) { FALSE })
setMethod("supported_constraints", "Elemental", function(object) { c(supported_constraints(ConicSolver()), "SOC") })

# TODO: Map of Elemental status to CVXR status.
setMethod("status_map", "Elemental", function(object, status) {
  if(status == 0)
    return(OPTIMAL)
  else
    stop("Unimplemented")
})

setMethod("import_solver", "Elemental", function(object) {
  stop("Elemental is currently unavailable in R. Please see http://libelemental.org/about/ for more information.")
})

setMethod("name", "Elemental", function(x) { ELEMENTAL_NAME })
setMethod("accepts", signature(object = "Elemental", problem = "Problem"), function(object, problem) {
  # Can Elemental solve the problem?
  # TODO: Check if the matrix is stuffed.
  if(!is_affine(problem@objective@args[[1]]))
    return(FALSE)
  for(constr in problem@constraints) {
    if(!class(constr) %in% supported_constraints(object))
      return(FALSE)
    for(arg in constr@args) {
      if(!is_affine(arg))
        return(FALSE)
    }
  }
  return(TRUE)
})

setMethod("apply", signature(object = "Elemental", problem = "Problem"), function(object, problem) {
  data <- list()
  objective <- canonical_form(problem@objective)[[1]]
  constraints <- lapply(problem@constraints, function(c) { canonical_form(c)[[2]] })
  constraints <- unlist(constraints, recursive = TRUE)
  data$objective <- objective
  data$constraints <- constraints
  inv_data <- list()
  inv_data[var_id(object)] <- id(variables(problem)[[1]])

  # Order and group constraints.
  eq_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "Zero" })]
  inv_data[eq_constr(object)] <- eq_constr
  leq_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "NonPos" })]
  inv_data[neq_constr(object)] <- leq_constr
  return(list(data, inv_data))
})

setMethod("solve_via_data", "Elemental", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  solver <- EL_OLD()
  prob_data <- list()
  prob_data[name(object)] <- ProblemData()
  solve(solver, data$objective, data$constraints, prob_data, warm_start, verbose, solver_opts)
})

GLPK <- setClass("GLPK", contains = "CVXOPT")
setMethod("mip_capable", "GLPK", function(object) { FALSE })
setMethod("supported_constraints", "GLPK", function(object) { supported_constraints(ConicSolver()) })

setMethod("name", "GLPK", function(x) { GLPK_NAME })
setMethod("import_solver", "GLPK", function(object) {
  requireNamespace("Rglpk", quietly = TRUE)
})

setMethod("accepts", signature(object = "GLPK", problem = "Problem"), function(object, problem) {
  # TODO: Check if the matrix is stuffed.
  if(!is_affine(problem@objective@args[[1]]))
    return(FALSE)
  for(constr in problem@constraints) {
    if(!class(constr) %in% supported_constraints(object))
      return(FALSE)
    for(arg in constr@args) {
      if(!is_affine(arg))
        return(FALSE)
    }
  }
  return(TRUE)
})

setMethod("invert", signature(object = "GLPK", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  status <- solution$status

  primal_vars <- NA
  dual_vars <- NA
  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution$value
    primal_vars <- list()
    primal_vars[inverse_data[var_id(object)]] <- solution$primal
  } else {
    if(status == INFEASIBLE)
      opt_val <- Inf
    else if(status == UNBOUNDED)
      opt_val <- -Inf
    else
      opt_val <- NA
  }

  return(Solution(status, opt_val, primal_vars, dual_vars, list()))
})

setMethod("solve_via_data", "GLPK", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  solver <- GLPK_OLD()
  prob_data <- list()
  prob_data[name(object)] <- ProblemData()
  solve(solver, data$objective, data$constraints, prob_data, warm_start, verbose, solver_opts)
})

GLPK_MI <- setClass("GLPK_MI", contains = "GLPK")
setMethod("mip_capable", "GLPK_MI", function(object) { TRUE })
setMethod("supported_constraints", "GLPK_MI", function(object) { supported_constraints(ConicSolver()) })
setMethod("name", "GLPK_MI", function(x) { GLPK_MI_NAME })
setMethod("solve_via_data", "GLPK_MI", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  solver <- GLPK_OLD()
  solver_opts[BOOL_IDX] <- data[BOOL_IDX]
  solver_opts[INT_IDX] <- data[INT_IDX]
  prob_data <- list()
  prob_data[name(object)] <- ProblemData()
  solve(solver, data$objective, data$constraints, prob_data, warm_start, verbose, solver_opts)
})

GUROBI <- setClass("GUROBI", contains = "ConicSolver")

# Solver capabilities.
setMethod("mip_capable", "GUROBI", function(object) { TRUE })
setMethod("supported_constraints", "GUROBI", function(object) { c(supported_constraints(ConicSolver()), "SOC") })

# Map of Gurobi status to CVXR status.
setMethod("status_map", "GUROBI", function(object, status) {
  if(status == 2)
    return(OPTIMAL)
  else if(status == 3)
    return(INFEASIBLE)
  else if(status == 5)
    return(UNBOUNDED)
  else if(status %in% c(4, 6, 7, 8, 10, 11, 12, 13))
    return(SOLVER_ERROR)
  else if(status == 9)   # TODO: Could be anything. Means time expired.
    return(OPTIMAL_INACCURATE)
  else
    stop("GUROBI status unrecognized: ", status)
})

setMethod("name", "GUROBI", function(x) { GUROBI_NAME })
setMethod("import_solver", "GUROBI", function(object) {
  requireNamespace("gurobi", quietly = TRUE)
})

setMethod("accepts", signature(object = "GUROBI", problem = "Problem"), function(object, problem) {
  # TODO: Check if the matrix is stuffed.
  if(!is_affine(problem@objective@args[[1]]))
    return(FALSE)
  for(constr in problem@constraints) {
    if(!class(constr) %in% supported_constraints(object))
      return(FALSE)
    for(arg in constr@args) {
      if(!is_affine(arg))
        return(FALSE)
    }
  }
  return(TRUE)
})

setMethod("apply", signature(object = "GUROBI", problem = "Problem"), function(object, problem) {
  data <- list()
  objective <- canonical_form(problem@objective)[[1]]
  constraints <- lapply(problem@constraints, function(c) { canonical_form(c)[[2]] })
  constraints <- unlist(constraints, recursive = TRUE)
  data$objective <- objective
  data$constraints <- constraints
  variables <- variables(problem)[[1]]
  data[BOOL_IDX] <- lapply(variables@boolean_idx, function(t) { t[1] })
  data[INT_IDX] <- lapply(variables@integer_idx, function(t) { t[1] })

  # Order and group constraints.
  inv_data <- list()
  inv_data[var_id(object)] <- id(variables(problem)[[1]])
  eq_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "Zero" })]
  inv_data[eq_constr(object)] <- eq_constr
  leq_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "NonPos" })]
  soc_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "SOC" })]
  inv_data[neq_constr(object)] <- c(leq_constr, soc_constr)
  inv_data$is_mip <- length(data[BOOL_IDX]) > 0 || length(data[INT_IDX]) > 0
  return(list(data, inv_data))
})

setMethod("invert", signature(object = "GUROBI", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  status <- solution$status

  primal_vars <- NA
  dual_vars <- NA
  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution$value
    primal_vars <- list()
    primal_vars[inverse_data[var_id(object)]] <- solution$primal
    if(!inverse_data$is_mip) {
      eq_dual <- get_dual_values(solution$eq_dual, extract_dual_value, inverse_data(eq_constr(object)))
      leq_dual <- get_dual_values(solution$ineq_dual, extract_dual_value, inverse_data(neq_constr(object)))
      eq_dual <- modifyList(eq_dual, leq_dual)
      dual_vars <- eq_dual
    }
  } else {
    if(status == INFEASIBLE)
      opt_val <- Inf
    else if(status == UNBOUNDED)
      opt_val <- -Inf
    else
      opt_val <- NA
  }

  return(Solution(status, opt_val, primal_vars, dual_vars, list()))
})

setMethod("solve_via_data", "GUROBI", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  solver <- GUROBI_OLD()
  solver_opts[BOOL_IDX] <- data[BOOL_IDX]
  solver_opts[INT_IDX] <- data[INT_IDX]
  prob_data <- list()
  prob_data[name(object)] <- ProblemData()
  solve(solver, data$objective, data$constraints, prob_data, warm_start, verbose, solver_opts)
})

LS <- setClass("LS", contains = "Solver")

# LS is incapable of solving any general cone program
# and must be invoked through a special path.
setMethod("lp_capable", "LS", function(object) { FALSE })
setMethod("socp_capable", "LS", function(object) { FALSE })
setMethod("psd_capable", "LS", function(object) { FALSE })
setMethod("exp_capable", "LS", function(object) { FALSE })
setMethod("mip_capable", "LS", function(object) { FALSE })

setMethod("import_solver", "LS", function(object) { })
setMethod("name", "LS", function(x) { LS_NAME })
setMethod("split_constr", "LS", function(object, constr_map) {
  return(list(constr_map[EQ_MAP], constr_map[LEQ_MAP], list()))
})

setMethod("suitable", signature(object = "LS", problem = "Problem"), function(object, problem) {
  allowedVariables <- c("Variable")
  # TODO: Handle affine objective.
  is_dcp(problem) &&
      is_quadratic(problem@objective@args[[1]]) &&
      !is_affine(problem@objective@args[[1]]) &&
      all(sapply(problem@constraints, function(c) { is(c, "Zero") })) &&
      all(sapply(variables(problem), function(v) { class(v) %in% allowedVariables })) &&
      all(sapply(variables(problem), function(v) { is.na(domain(v)) || length(domain(v)) == 0 }))   # No implicit variable domains
  # TODO: Domains are not implemented yet.
})

setMethod("validate_solver", signature(object = "LS", problem = "Problem"), function(object, problem) {
  if(!suitable(object, problem))
    stop("The solver ", name(object), " cannot solve the problem.")
})

setMethod("get_sym_data", "LS", function(solver, objective, constraints, cached_data = NA) {
  FakeSymData <- setClass("FakeSymData", representation(constr_map = "list", vars_ = "list", var_offsets = "numeric", var_sizes = "list", x_length = "numeric"),
                                         prototype(constr_map = list(), vars_ = list(), var_offsets = NA_integer_, var_sizes = list(), x_length = NA_real_))

  setMethod("initialize", "FakeSymData", function(.Object) {
    .Object@constr_map <- list()
    .Object@constr_map[EQ_MAP] <- constraints
    vars_ <- variables(objective)
    for(c in constraints)
      vars_ <- c(vars_, variables(c))
    vars_ <- unique(vars_)
    .Object@vars_ <- vars_
    offsets <- get_var_offsets(solver, vars_)
    .Object@var_offsets <- offsets[[1]]
    .Object@var_sizes <- offsets[[2]]
    .Object@x_length <- offsets[[2]]
    return(.Object)
  })

  setMethod("get_var_offsets", "FakeSymData", function(object, variables) {
    var_offsets <- list()
    var_sizes <- list()
    vert_offset <- 0
    for(x in variables) {
      x_id <- as.character(id(x))
      var_sizes[x_id] <- size(x)
      var_offsets[x_id] <- vert_offset
      vert_offset <- vert_offset + size(x)[1]*size(x)[2]
    }
    return(list(var_offsets, var_sizes, vert_offset))
  })

  return(FakeSymData(objective, constraints))
})

setMethod("solve", "LS", function(solver, objective, constraints, cached_data, warm_start, verbose, solver_opts) {
  sym_data <- get_sym_data(solver, objective, constraints)

  id_map <- sym_data@var_offsets
  N <- sym_data@x_length
  extractor <- CoeffExtractor(id_map, N)

  # Extract the coefficients.
  coeffs <- get_coeffs(extractor, objective@args[[1]])
  Ps <- coeffs[[1]]
  Q <- coeffs[[2]]
  R <- coeffs[[3]]

  P <- Ps[1]
  q <- as.vector(Q)
  r <- R[1]

  # Forming the KKT system.
  if(length(constraints) > 0) {
    Cs <- lapply(constraints, function(c) { get_coeffs(constraints, c@.expr)[-1] })
    As <- do.call("rbind", lapply(Cs, function(C) { C[1] }))
    bs <- sapply(Cs, function(C) { C[2] })
    lhs <- rbind(cbind(2*P, t(As)), cbind(As, matrix(0, nrow = ncol(As), ncol = nrow(As))))
    lhs <- Matrix(lhs, sparse = TRUE)
    rhs <- c(-q, -bs)
  } else {
    lhs <- 2*P
    rhs <- -q
  }

  # TODO: Actually solve the KKT system.
  tryCatch({
    sol <- SLA.spsolve(lhs, rhs)
    x <- sol[1:N]
    nu <- sol[(N+1):length(sol)]
    p_star <- t(x) %*% (P %*% x + q) + r
  }, error = function(e) {
    if(e == SLA.MatrixRankWarning) {
      x <- NA
      nu <- NA
      p_star <- NA
    }
  })

  result_dict <- list()
  result_dict[PRIMAL] <- x
  result_dict[DUAL] <- nu
  result_dict[VALUE] <- p_star
  return(format_results(result_dict, NA, cached_data))
})

setMethod("format_results", "LS", function(object, result_dict, data, cached_data) {
  new_results <- result_dict
  if(is.na(result_dict[PRIMAL]) || is.null(result_dict[PRIMAL]))
    new_results[STATUS] <- INFEASIBLE
  else
    new_results[STATUS] <- OPTIMAL
  return(new_results)
})

vectorized_lower_tri_to_mat <- function(v, dim) {
  # Return the symmetric 2D array defined by taking "v" to specify its lower triangular entries.
  rows <- c()
  cols <- c()
  vals <- c()
  running_idx <- 1
  for(j in 1:dim) {
    rows <- c(rows, j + 0:(dim - j))
    cols <- c(cols, rep(j, dim - j + 1))
    vals <- c(vals, v[running_idx:(running_idx + dim - j + 1)])
    running_idx <- running_idx + dim - j + 1
  }
  A <- sparseMatrix(i = rows, j = cols, x = vals, dims = c(dim, dim))
  d <- diag(diag(A))
  A <- A + t(A) - d
  return(A)
}

psd_coeff_offset <- function(problem, c) {
  # Returns an array "G" and vector "h" such that the given constraint is equivalent to "G * z <=_{PSD} h".
  extractor <- CoeffExtractor(InverseData(problem))
  Ab_vec <- affine(extractor, c@expr)
  G <- -Ab_vec[[1]]
  h <- Ab_vec[[2]]
  dim <- nrow(c@expr)
  return(list(G, h, dim))
}

MOSEK <- setClass("MOSEK", contains = "ConicSolver")

setMethod("mip_capable", "MOSEK", function(object) { TRUE })
setMethod("supported_constraints", "MOSEK", function(object) { c(supported_constraints(ConicSolver()), "SOC", "PSD") })
setMethod("exp_cone_order", "MOSEK", function(object) { c(2, 1, 0) })

setMethod("import_solver", "MOSEK", function(object) {
  requireNamespace("Rmosek", quietly = TRUE)
  # TODO: Add exponential cone support.
})

setMethod("name", "MOSEK", function(x) { MOSEK_NAME })
setMethod("accepts", signature(object = "MOSEK", problem = "Problem"), function(object, problem) {
  # TODO: Check if the matrix is stuffed.
  import_solver(object)
  if(!is_affine(problem@objective@args[[1]]))
    return(FALSE)
  for(constr in problem@constraints) {
    if(!class(constr) %in% supported_constraints(object))
      return(FALSE)
    for(arg in constr@args) {
      if(!is_affine(arg))
        return(FALSE)
    }
  }
  return(TRUE)
})

setMethod("block_format", "MOSEK", function(object, problem, constraints, exp_cone_order = NA) {
  if(length(constraints) == 0 || is.null(constraints) || is.na(constraints))
    return(list(NA, NA))
  matrices <- list()
  offsets <- list()
  lengths <- c()
  ids <- c()

  for(con in constraints) {
    coeff_offs <- format_constr(object, problem, con, exp_cone_order)
    coeff <- coeff_offs[[1]]
    offset <- coeff_offs[[2]]
    matrices <- c(matrices, coeff)
    offsets <- c(offsets, offset)
    lengths <- c(size(offset))
    ids <- c(ids, id(con))
  }
  coeff <- Matrix(do.bind("rbind", matrices), sparse = TRUE)
  offset <- do.bind("cbind", offsets)
  return(list(coeff, offset, lengths, ids))
})

setMethod("apply", signature(object = "MOSEK", problem = "Problem"), function(object, problem) {
  data <- list()
  inv_data <- list(suc_slacks = list(), y_slacks = list(), snx_slacks = list(), psd_dims = list())
  inv_data[var_id(object)] <- id(variables(problem)[[1]])

  # Get integrality constraint information.
  var <- variables(problem)[[1]]
  data[BOOL_IDX] <- sapply(var@boolean_idx, function(t) { as.integer(t[1]) })
  data[INT_IDX] <- sapply(var@integer_idx, function(t) { as.integer(t[1]) })
  inv_data$integer_variables <- length(data[BOOL_IDX]) + length(data[INT_IDX]) > 0

  # Parse the coefficient vector from the objective.
  coeff_offs <- get_coeff_offset(object, problem@objective@args[[1]])
  c <- coeff_offs[[1]]
  constant <- coeff_offs[[2]]
  data[C_MAP] <- as.vector(c)
  inv_data$n0 <- length(data[C_KEY])
  data[OBJ_OFFSET] <- constant[1]
  data[DIMS] <- list()
  data[DIMS][SOC_DIM] <- list()
  data[DIMS][EXP_DIM] <- list()
  data[DIMS][PSD_DIM] <- list()
  data[DIMS][LEQ_DIM] <- 0
  data[DIMS][EQ_DIM] <- 0
  inv_data[OBJ_OFFSET] <- constant[1]
  Gs <- list()
  hs <- list()

  # Linear inequalities.
  leq_constr <- problem@constraints[sapply(problem@constraints, function(ci) { class(ci) == "NonPos" })]
  if(length(leq_constr) > 0) {
    blform <- block_format(object, problem, leq_constr)   # G, h : G*z <= h.
    G <- blform[[1]]
    h <- blform[[2]]
    lengths <- blform[[3]]
    ids <- blform[[4]]
    inv_data$suc_slacks <- c(inv_data$suc_slacks, lapply(1:length(lengths), function(k) { c(ids[k], lengths[k]) }))
    data[DIMS][LEQ_DIM] <- sum(lengths)
    Gs <- c(Gs, G)
    hs <- c(hs, h)
  }

  # Linear equations.
  eq_constr <- problem@constraints[sapply(problem@constraints, function(ci) { class(ci) == "Zero" })]
  if(length(leq_constr) > 0) {
    blform <- block_format(object, problem, eq_constr)   # G, h : G*z == h.
    G <- blform[[1]]
    h <- blform[[2]]
    lengths <- blform[[3]]
    ids <- blform[[4]]
    inv_data$y_slacks <- c(inv_data$y_slacks, lapply(1:length(lengths), function(k) { c(ids[k], lengths[k]) }))
    data[DIMS][EQ_DIM] <- sum(lengths)
    Gs <- c(Gs, G)
    hs <- c(hs, h)
  }

  # Second order cone.
  soc_constr <- problem@constraints[sapply(problem@constraints, function(ci) { class(ci) == "SOC" })]
  data[DIMS][SOC_DIM] <- list()
  for(ci in soc_constr)
    data[DIMS][SOC_DIM] <- c(data[DIMS][SOC_DIM], cone_sizes(ci))
  if(length(soc_constr) > 0) {
    blform <- block_format(object, problem, soc_constr)   # G*z <=_{soc} h.
    G <- blform[[1]]
    h <- blform[[2]]
    lengths <- blform[[3]]
    ids <- blform[[4]]
    inv_data$snx_slacks <- c(inv_data$snx_slacks, lapply(1:length(lengths), function(k) { c(ids[k], lengths[k]) }))
    Gs <- c(Gs, G)
    hs <- c(hs, h)
  }

  # Exponential cone.
  exp_constr <- problem@constraints[sapply(problem@constraints, function(ci) { class(ci) == "ExpCone" })]
  if(length(exp_constr) > 0) {
    # G*z <=_{EXP} h.
    blform <- block_format(object, problem, exp_constr, exp_cone_order(object))
    G <- blform[[1]]
    h <- blform[[2]]
    lengths <- blform[[3]]
    ids <- blform[[4]]
    data[DIMS][EXP_DIM] <- lengths
    Gs <- c(Gs, G)
    hs <- c(hs, h)
  }

  # PSD constraints.
  psd_constr <- problem@constraints[sapply(problem@constraints, function(ci) { class(ci) == "PSD" })]
  if(length(psd_constr) > 0) {
    data[DIMS][PSD_DIM] <- list()
    for(c in psd_constr) {
      coeff_offs <- psd_coeff_offset(problem, c)
      G_vec <- coeff_offs[[1]]
      h_vec <- coeff_offs[[2]]
      dim <- coeff_offs[[3]]
      inv_data$psd_dims <- c(inv_data$psd_dims, c(id(c), dim))
      data[DIMS][PSD_DIM] <- c(data[DIMS][PSD_DIM], dim)
      Gs <- c(Gs, G_vec)
      hs <- c(hs, h_vec)
    }
  }

  data[G_KEY] <- Matrix(do.call("rbind", Gs), sparse = TRUE)
  data[H_KEY] <- Matrix(do.call("cbind", hs), sparse = TRUE)
  inv_data$is_LP <- (length(psd_constr) + length(exp_constr) + length(soc_constr)) == 0
  return(list(data, inv_data))
})

# TODO: Finish MOSEK class implementation.
setMethd("solve_via_data", "MOSEK", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  requireNamespace("Rmosek", quietly = TRUE)
  env <- Rmosek::Env()
  task <- env.Task(0,0)

  # TODO: Handle logging for verbose.

  # Parse all user-specified parameters (override default logging
  # parameters if applicable).
  kwargs <- sort(names(solver_opts))
  save_file <- NA_character_
  if("mosek_params" %in% kwargs) {
    MOSEK._handle_mosek_params(task, solver_opts$mosek_params)
    kwargs$mosek_params <- NULL
  }
  if("save_file" %in% kwargs) {
    save_file <- solver_opts$save_file
    kwargs$save_file <- NULL
  }
  if("bfs" %in% kwargs)
    kwargs$bfs <- NULL
  if(!(is.null(kwargs) || length(kwargs) == 0))
    stop("Invalid keyword-argument ", kwargs[[1]])

  # Check if the CVXR standard form has zero variables. If so,
  # return a trivial solution. This is necessary because MOSEK
  # will crash if handed a problem with zero variables.
  if(length(data[C_KEY]) == 0) {
    res <- list()
    res[STATUS] <- OPTIMAL
    res[PRIMAL] <- list()
    res[VALUE] <- data[OFFSET]
    res[EQ_DUAL] <- list()
    res[INEQ_DUAL] <- list()
    return(res)
  }

  # The following lines recover problem parameters, and define helper constants.
  #
  #   The problem's objective is "min c.T * z".
  #   The problem's constraint set is "G * z <=_K h."
  #   The rows in (G, h) are formatted in order of
  #       (1) linear inequalities,
  #       (2) linear equations,
  #       (3) soc constraints,
  #       (4) exponential cone constraints,
  #       (5) vectorized linear matrix inequalities.
  #   The parameter "dims" indicates the exact
  #   dimensions of each of these cones.
  #
  #   MOSEK's standard form requires that we replace generalized
  #   inequalities with slack variables and linear equations.
  #   The parameter "n" is the size of the column-vector variable
  #   after adding slacks for SOC and EXP constraints. To be
  #   consistent with MOSEK documentation, subsequent comments
  #   refer to this variable as "x".

  c <- data[C_KEY]
  G <- data[G_KEY]
  h <- data[H_KEY]
  n0 <- length(c)
  n <- n0 + sum(dims[SOC_DIM]) + sum(dims[EXP_DIM])
  psd_total_dims <- sum(dims[PSD_DIM]^2)
  m <- length(h)
  num_bool <- length(data[BOOL_IDX])
  num_int <- length(data[INT_IDX])

  # Define variables, cone constraints, and integrality constraints.
  #
  #   The variable "x" is a length-n block vector, with
  #       Block 1: "z" from "G * z <=_K h",
  #       Block 2: slacks for SOC constraints, and
  #       Block 3: slacks for EXP cone constraints.
  #
  #   Once we declare x in the MOSEK model, we add the necessary
  #   conic constraints for slack variables (Blocks 2 and 3).
  #   The last step is to add integrality constraints.
  #
  #   Note that the API call for PSD variables contains the word "bar".
  #   MOSEK documentation consistently uses "bar" as a sort of flag,
  #   indicating that a function deals with PSD variables.

  task.appendvars(n)
  task.putvarboundlist(1:n, rep(mosek.boundkey.fr, n), matrix(0, nrow = n, ncol = 1), matrix(0, nrow = n, ncol = 1))
  if(psd_total_dims > 0)
    task.appendbarvars(dims[PSD_DIM])
  running_idx <- n0
  for(size_cone in dims[SOC_DIM]) {
    task.appendcone("quad", 0, running_idx:(running_idx + size_cone))
    running_idx <- running_idx + size_cone
  }
  for(k in 1:floor(sum(dims[EXP_DIM])/3)) {
    task.appendcone("pexp", 0, running_idx:(running_idx + 3))
    running_idx <- running_idx + 3
  }
  if(num_bools + num_int > 0) {
    task.putvartypelist(data[BOOL_IDX], rep("integer", num_bool))
    task.putvarboundlist(data[BOOL_IDX], rep(boundkey, num_bool), rep(0, num_bool), rep(1, num_bool))
    task.putvartypelist(data[INT_IDX], rep("integer", num_int))
  }

  # Define linear inequality and equality constraints.
  #
  #   Mosek will see a total of m linear expressions, which must
  #   define linear inequalities and equalities. The variable x
  #   contributes to these linear expressions by standard
  #   matrix-vector multiplication; the matrix in question is
  #   referred to as "A" in the mosek documentation. The PSD
  #   variables have a different means of contributing to the
  #   linear expressions. Specifically, a PSD variable Xj contributes
  #   "+tr( \bar{A}_{ij} * Xj )" to the i-th linear expression,
  #   where \bar{A}_{ij} is specified by a call to putbaraij.
  #
  #   The following code has three phases.
  #       (1) Build the matrix A.
  #       (2) Specify the \bar{A}_{ij} for PSD variables.
  #       (3) Specify the RHS of the m linear (in)equalities.
  #
  #   Remark : The parameter G gives every row in the first
  #   n0 columns of A. The remaining columns of A are for SOC
  #   and EXP slack variables. We can actually account for all
  #   of these slack variables at once by specifying a giant
  #   identity matrix in the appropriate position in A.

  task.appendcons(m)
  G_sparse <- as(as.matrix(G), "sparseMatrix")
  G_sum <- summary(G_sparse)
  rows <- G_sum$i
  cols <- G_sum$j
  vals <- G_sum$x
  task.putaijlist(as.list(row), as.list(col), as.list(vals))
  total_soc_exp_slacks <- sum(dims[SOC_DIM]) + sum(dims[EXP_DIM])
  if(total_soc_exp_slacks > 0) {
    i <- dims[LEQ_DIM] + dims[EQ_DIM]   # Constraint index in (1, ..., m)
    j <- length(c)   # Index of the first slack variable in the block vector "x".
    rows <- as.list(i:(i + total_soc_exp_slacks))
    cols <- as.list(j:(j + total_soc_exp_slacks))
    task.putaijlist(rows, cols, rep(1, total_soc_exp_slacks))
  }

  # Constraint index: start of LMIs.
  i <- dims[LEQ_DIM] + dims[EQ_DIM] + total_soc_exp_slacks
  for(j in 1:dims[PSD_DIM]) {   # SDP slack variable "Xj"
    dim <- dims[PSD_DIM][j]
    for(row_idx in 1:dim) {
      for(col_idx in 1:dim) {
        val <- ifelse(row_idx == col_idx, 1, 0.5)
        row <- max(row_idx, col_idx)
        col <- min(row_idx, col_idx)
        mat <- task.appendsparsesymmat(dim, list(row), list(col), list(val))
        task.putbaraij(i, j, list(mat), list(1.0))
        i <- i + 1
      }
    }
  }

  num_eq <- length(h) - dims[LEQ_DIM]
  type_constraint <- rep(mosek.boundkey.up, dims[LEQ_DIM]) + rep(mosek.boundkey.fx, num_eq)
  task.putconboundlist(1:m, type_constraint, h, h)

  # Define the objective and optimize the MOSEK task.
  task.putclist(1:length(c), c)
  task.putobjsense(mosek.objsense.minimize)
  if(!is.na(save_file))
    task.writedata(save_file)
  task.optimize

  if(verbose)
    task.solutionsummary(mosek.streamtype.msg)

  return(list(env = env, task = task, solver_options = solver_opts))
})

setMethod("invert", "MOSEK", function(object, results, inverse_data) {
  requireNamespace("Rmosek", quietly = TRUE)

  has_attr <- !is.null(mosek.solsta$near_optimal)
  status_map <- function(status) {
    if(status %in% c("optimal", "integer_optimal"))
      return(OPTIMAL)
    else if(status %in% c("prim_feas", "near_optimal", "near_integer_optimal"))
      return(OPTIMAL_INACCURATE)
    else if(status == "prim_infeas_cer") {
      if(has_attr)
        return(INFEASIBLE_INACCURATE)
      else
        return(INFEASIBLE)
    } else if(status == "dual_infeas_cer") {
      if(has_attr)
        return(UNBOUNDED_INACCURATE)
      else
        return(UNBOUNDED)
    } else
      return(SOLVER_ERROR)
  }

  env <- results$env
  task <- results$task
  solver_opts <- results$solver_options

  if(inverse_data$integer_variables)
    sol <- mosek.soltype.itg
  else if(!is.null(solver_opts$bfs) && solver_opts$bfs && inverse_data$is_LP)
    sol <- mosek.soltype.bas   # The basic feasible solution.
  else
    sol <- mosek.soltype.itr   # The solution found via interior point method.

  problem_status <- task.getprosta(sol)
  solution_status <- task.getsolsta(sol)

  status <- status_map(solution_status)

  # For integer problems, problem status determines infeasibility (no solution).
  if(sol == mosek.soltype.itg && problem_status == mosek.prosta.prim_infeas)
    status <- INFEASIBLE

  if(status %in% SOLUTION_PRESENT) {
    # Get objective value.
    opt_val <- task.getprimalobj(sol) + inverse_data[OBJ_OFFSET]
    # Recover the CVXR standard form primal variable.
    z <- rep(0, inverse_data$n0)
    task.getxxslice(sol, 0, length(z), z)
    primal_vars <- list()
    primal_vars[inverse_data[var_id(object)]] <- z
    # Recover the CVXR standard form dual variables.
    if(sol == mosek.soltype.itg)
      dual_vars <- NA
    else
      dual_vars <- MOSEK.recover_dual_variables(task, sol, inverse_data)
  } else {
    if(status == INFEASIBLE)
      opt_val <- Inf
    else if(status == UNBOUNDED)
      opt_val <- -Inf
    else
      opt_val <- NA
    primal_vars <- NA
    dual_vars <- NA
  }

  # Store computation time.
  attr <- list()
  attr[SOLVE_TIME] <- task.getdouinf(mosek.dinfitem.optimizer_time)

  # Delete the MOSEK Task and Environment
  task.__exit__(NA, NA, NA)
  env.__exit__(NA, NA, NA)

  return(Solution(status, opt_val, primal_vars, dual_vars, attr))
})

MOSEK.recover_dual_variables <- function(task, sol, inverse_data) {
  dual_vars <- list()

  # Dua variables for the inequality constraints.
  suc_len <- sum(sapply(inverse_data$suc_slacks, function(val) { val[[2]] }))
  if(suc_len > 0) {
    suc <- rep(0, suc_len)
    task.getsucslice(sol, 0, suc_len, suc)
    dual_vars <- modifyList(dual_vars, MOSEK.parse_dual_vars(suc, inverse_data$suc_slacks))
  }

  # Dual variables for the original equality constraints.
  y_len <- sum(sapply(inverse_data$y_slacks, function(val) { val[[2]] }))
  if(y_len > 0) {
    y <- rep(0, y_len)
    task.getyslice(sol, suc_len, suc_len + y_len, y)
    dual_vars <- modifyList(dual_vars, MOSEK.parse_dual_vars(y, inverse_data$y_slacks))
  }

  # Dual variables for SOC and EXP constraints.
  snx_len <- sum(sapply(inverse_data$snx_slacks, function(val) { val[[2]] }))
  if(snx_len > 0) {
    snx <- matrix(0, nrow = snx_len, ncol = 1)
    task.getsnxslice(sol, inverse_data$n0, inverse_data$n0 + snx_len, snx)
    dual_vars <- modifyList(dual_vars, MOSEK.parse_dual_vars(snx, inverse_data$snx_slacks))
  }

  # Dual variables for PSD constraints.
  for(j in 1:length(inverse_data$psd_dims)) {
    id <- inverse_data$psd_dims[j][[1]]
    dim <- inverse_data$psd_dims[j][[2]]
    sj <- rep(0, dim*floor((dim + 1)/2))
    task.getbars(sol, j, sj)
    dual_vars[id] <- vectorized_lower_tri_to_mat(sj, dim)
  }
  return(dual_vars)
}

MOSEK.parse_dual_vars <- function(dual_var, constr_id_to_constr_dim) {
  dual_vars <- list()
  running_idx <- 1
  for(val in constr_id_to_constr_dim) {
    id <- val[[1]]
    dim <- val[[2]]
    if(dim == 1)
      dual_vars[id] <- dual_vars[running_idx]   # a scalar.
    else
      dual_vars[id] <- as.matrix(dual_vars[running_idx:(running_idx + dim)])
    running_idx <- running_idx + dim
  }
  return(dual_vars)
}

MOSEK._handle_mosek_params <- function(task, params) {
  if(is.na(params))
    return()

  requireNamespace("Rmosek", quietly = TRUE)

  handle_str_param <- function(param, value) {
    if(startsWith(param, "MSK_DPAR_"))
      task.putnadourparam(param, value)
    else if(startsWith(param, "MSK_IPAR_"))
      task.putnaintparam(param, value)
    else if(startsWith(param, "MSK_SPAR_"))
      task.putnastrparam(param, value)
    else
      stop("Invalid MOSEK parameter ", param)
  }

  handle_enum_param <- function(param, value) {
    if(is(param, "dparam"))
      task.putdouparam(param, value)
    else if(is(param, "iparam"))
      task.putintparam(param, value)
    else if(is(param, "sparam"))
      task.putstrparam(param, value)
    else
      stop("Invalid MOSEK parameter ", param)
  }

  for(p in params) {
    param <- p[[1]]
    value <- p[[2]]
    if(is(param, "character"))
      handle_str_param(param, value)
    else
      handle_enum_param(param, value)
  }
}

# Utility method for formatting a ConeDims instance into a dictionary
# that can be supplied to SCS.
dims_to_solver_dict <- function(cone_dims) {
  cones <- list(f = as.integer(cone_dims@zero),
                l = as.integer(cone_dims@nonpos),
                q = sapply(cone_dims@soc, as.integer),
                ep = as.integer(cone_dims@exp),
                s = sapply(cone_dims@psd, as.integer))
  return(cones)
}

# Utility methods for special handling of semidefinite constraints.
scaled_lower_tri <- function(matrix) {
  # Returns an expression representing the lower triangular entries.
  # Scales the strictly lower triangular entries by sqrt(2), as
  # required by SCS.
  rows <- cols <- nrow(matrix)
  entries <- rows * floor((cols + 1)/2)
  val_arr <- c()
  row_arr <- c()
  col_arr <- c()
  count <- 0

  for(j in 1:cols) {
    for(i in 1:rows) {
      if(j <= i) {
        # Index in the original matrix.
        col_arr <- c(col_arr, (j-1)*rows + i)
        # Index in the extracted vector.
        row_arr <- c(row_arr, count + 1)
        if(j == i)
          val_arr <- c(val_arr, 1.0)
        else
          val_arr <- c(val_arr, sqrt(2))
        count <- count + 1
      }
    }
  }

  shape <- c(entries, rows*cols)
  coeff <- Constant(sparseMatrix(i = row_arr, j = col_arr, x = val_arr), shape)
  vectorized_matrix <- reshape(matrix, c(rows*cols, 1))
  return(coeff %*% vectorized_matrix)
}

tri_to_full <- function(lower_tri, n) {
  # Expands n*floor((n+1)/2) lower triangular to full matrix.
  # Scales off-diagonal by 1/sqrt(2), as per the SCS specification.
  full <- matrix(0, nrow = n, ncol = n)
  for(col in 1:n) {
    for(row in col:n) {
      idx <- row - col + n*floor((n+1)/2) - (n-col+1)*floor((n-col+2)/2) + 1
      if(row != col) {
        full[row, col] <- lower_tri[idx]/sqrt(2)
        full[col, row] <- lower_tri[idx]/sqrt(2)
      } else
        full[row, col] <- lower_tri[idx]
    }
  }
  matrix(full, nrow = n*n, byrow = FALSE)
}

SCS <- setClass("SCS", contains = "ConicSolver")

# Solver capabilities.
setMethod("mip_capable", "SCS", function(object) { FALSE })
setMethod("supported_constraints", "SCS", function(object) { c(supported_constraints(ConicSolver()), "SOC", "ExpCone", "PSD") })

requires_constr.SCS <- function(object) { TRUE }

# Map of SCS status to CVXR status.
setMethod("status_map", "SCS", function(object, status) {
  if(status == "Solved")
    return(OPTIMAL)
  else if(status == "Solved/Inaccurate")
    return(OPTIMAL_INACCURATE)
  else if(status == "Unbounded")
    return(UNBOUNDED)
  else if(status == "Unbounded/Inaccurate")
    return(UNBOUNDED_INACCURATE)
  else if(status == "Infeasible")
    return(INFEASIBLE)
  else if(status == "Infeasible/Inaccurate")
    return(INFEASIBLE_INACCURATE)
  else if(status %in% c("Failure", "Indeterminate", "Interrupted"))
    return(SOLVER_ERROR)
  else
    stop("SCS status unrecognized: ", status)
})

# Order of exponential cone arguments for solver.
exp_cone_order.SCS <- function(object) { c(0, 1, 2) }

setMethod("name", "SCS", function(object) { SCS_NAME })
setMethod("import_solver", "SCS", function(object) {
  requireNamespace("scs", quietly = TRUE)
})

setMethod("format_constr", "SCS", function(object, problem, constr, exp_cone_order) {
  # Extract coefficient and offset vector from constraint.
  # Special cases PSD constraints, as SCS expects constraints to be
  # imposed on solely the lower triangular part of the variable matrix.
  # Moreover, it requires the off-diagonal coefficients to be scaled by
  # sqrt(2).
  if(is(constr, "PSD")) {
    expr <- constr@expr
    triangularized_expr <- scaled_lower_tri(expr)
    extractor <- CoeffExtractor(InverseData(problem))
    Ab <- affine(extractor, triangularized_expr)
    A_prime <- Ab[[1]]
    b_prim <- Ab[[2]]

    # SCS requests constraints to be formatted as Ax + s = b,
    # where s is constrained to reside in some cone. Here, however,
    # we are formatting the constraint as A"x + b" = -Ax + b; h ence,
    # A = -A", b = b".
    return(list(-1*A_prime, b_prime))
  } else
    callNextMethod(object, problem, constr, exp_cone_order)
})

setMethod("apply", signature(object = "SCS", problem = "Problem"), function(object, problem) {
  # Returns a new problem and data for inverting the new solution.
  data <- list()
  inv_data <- list()
  inv_data[var_id(object)] <- id(variables(problem)[[1]])

  # Parse the coefficient vector from the objective.
  offsets <- get_coeff_offset(object, problem@objective@args[[1]])
  data[C_KEY] <- offsets[[1]]
  data[OFFSET] <- offsets[[2]]
  data[C_KEY] <- as.vector(data[C_KEY])
  inv_data[OFFSET] <- data[OFFSET][1]

  # Order and group nonlinear constraints.
  constr_map <- group_constraints(problem@constraints)
  data[dims(ConicSolver())] <- ConeDims(constr_map)
  inv_data[dims(ConicSolver())] <- data[dims(ConicSolver())]

  # SCS requires constraints to be specified in the following order:
  # 1) Zero cone.
  # 2) Non-negative orthant.
  # 3) SOC.
  # 4) PSD.
  # 5) Exponential.
  zero_constr <- constr_map$Zero
  neq_constr <- c(constr_map$NonPos, constr_map$SOC, constr_map$PSD, constr_map$ExpCone)
  inv_data[eq_constr(object)] <- zero_constr
  inv_data[neq_constr(object)] <- neq_constr

  # Obtain A, b such that Ax + s = b, s \in cones.
  # Note that SCS mandates that the cones MUST be ordered with
  # zero cones first, then non-negative orthant, then SOC, then
  # PSD, then exponential.
  offsets <- group_coeff_offset(object, problem, c(zero_constr, neq_constr), exp_cone_order(object))
  data[A_KEY] <- offsets[[1]]
  data[B_KEY] <- offsets[[2]]
  return(list(data, inv_data))
})

setMethod("extract_dual_value", "SCS", function(object, result_vec, offset, constraint) {
  # Extracts the dual value for constraint starting at offset.
  # Special cases PSD constraints, as per the SCS specification.
  if(is(constraint, "PSD")) {
    dim <- nrow(constraint)
    lower_tri_dim <- dim*floor((dim+1)/2)
    new_offset <- offset + lower_tri_dim
    lower_tri <- result_vec[offset:new_offset]
    full <- tri_to_full(lower_tri, dim)
    return(list(full, new_offset))
  } else
    return(extract_dual_value(result_vec, offset, constraint))
})

setMethod("invert", signature(object = "SCS", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  # Returns the solution to the original problem given the inverse_data.
  status <- status_map(object, solution$info$status)

  attr <- list()
  attr[SOLVE_TIME] <- solution$info$solveTime
  attr[SETUP_TIME] <- solution$info$setupTime
  attr[NUM_ITERS] <- solution$info$iter

  if(status %in% SOLUTION_PRESENT) {
    primal_val <- solution$info$pobj
    opt_val <- primal_val + inverse_data[OFFSET]
    primal_vars <- list()
    primal_vars[inverse_data[var_id(object)]] <- as.matrix(solution$x)

    eq_dual_vars <- get_dual_values(as.matrix(solution$y[1:inverse_data[dims(ConicSolver())]@zero]),
      extract_dual_value(object), inverse_data[eq_constr(object)])

    ineq_dual_vars <- get_dual_values(as.matrix(solution$y[inverse_data[dims(ConicSolver())]@zero:length(solution$y)]),
                                      extract_dual_value(object), inverse_data[neq_constr(object)])

    dual_vars <- list()
    dual_vars <- modifyList(dual_vars, eq_dual_vars)
    dual_vars <- modifyList(dual_vars, ineq_dual_vars)
    return(Solution(status, opt_val, primal_vars, dual_vars, attr))
  } else
    return(failure_solution(status))
})

setMethod("solve_via_data", "SCS", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  # Returns the result of the call to the solver.
  requireNamespace("scs", quietly = TRUE)
  args <- list(A = data[A_KEY], b = data[B_KEY], c = data[C_KEY])
  if(warm_start && !is.na(solver_cache) && name(object) %in% names(solver_cache)) {
    args$x <- solver_cache[name(object)]$x
    args$y <- solver_cache[name(object)]$y
    args$s <- solver_cache[name(object)]$s
  }
  cones <- dims_to_solver_dict(data[dims(ConicSolver())])
  results <- scs::solve(args, cones, verbose = verbose, solver_opts)
  if(!is.na(solver_cache))
    solver_cache[name(object)] <- results
  return(results)
})

SuperSCS <- setClass("SuperSCS", contains = "SCS")
default_settings.SuperSCS <- function(object) {
  list(use_indirect = FALSE, eps = 1e-8, max_iters = 10000)
}

setMethod("name", "SuperSCS", function(x) { SUPER_SCS_NAME })
setMethod("import_solver", "SuperSCS", function(object) {
  stop("Unimplemented: SuperSCS is currently unavailable in R.")
})

setMethod("solve_via_data", "SuperSCS", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  args <- list(A = data[A_KEY], b = data[B_KEY], c = data[C_KEY])
  if(warm_start && !is.na(solver_cache) && name(object) %in% names(solver_cache)) {
    args$x <- solver_cache[name(object)]$x
    args$y <- solver_cache[name(object)]$y
    args$s <- solver_cache[name(object)]$s
  }
  cones <- dims_to_solver_dict(data[dims(ConicSolver())])

  # Settings.
  user_opts <- names(solver_opts)
  for(k in names(default_settings(SuperSCS))) {
    if(!k %in% user_opts)
      solver_opts[k] <- default_settings(SuperSCS)[k]
  }
  results <- SuperSCS::solve(args, cones, verbose = verbose, solver_opts)
  if(!is.na(solver_cache))
    solver_cache[name(object)] <- results
  return(results)
})

XPRESS <- setClass("XPRESS", contains = "ConicSolver")

# Solver capabilities.
setMethod("mip_capable", "XPRESS", function(object) { TRUE })
setMethod("supported_constraints", "XPRESS", function(object) { c(supported_constraints(ConicSolver()), "SOC") })

# Map of XPRESS status to CVXR status.
setMethod("status_map", "XPRESS", function(object, status) {
  if(status == 2)
    return(OPTIMAL)
  else if(status == 3)
    return(INFEASIBLE)
  else if(status == 5)
    return(UNBOUNDED)
  else if(status %in% c(4, 6, 7, 8, 10, 11, 12, 13))
    return(SOLVER_ERROR)
  else if(status == 9)   # TODO: Could be anything. Means time expired.
    return(OPTIMAL_INACCURATE)
  else
    stop("XPRESS status unrecognized: ", status)
})

setMethod("name", "XPRESS", function(x) { XPRESS_NAME })
setMethod("import_solver", "XPRESS", function(object) {
    stop("Unimplemented: XPRESS solver unavailable in R.")
})

setMethod("accepts", signature(object = "XPRESS", problem = "Problem"), function(object, problem) {
  # TODO: Check if the matrix is stuffed.
  if(!is_affine(problem@objective@args[[1]]))
    return(FALSE)
  for(constr in problem@constraints) {
    if(!class(constr) %in% supported_constraints(object))
      return(FALSE)
    for(arg in constr@args) {
      if(!is_affine(arg))
        return(FALSE)
    }
  }
  return(TRUE)
})

setMethod("apply", signature(object = "XPRESS", problem = "Problem"), function(object, problem) {
  data <- list()
  objective <- canonical_form(problem@objective)[[1]]
  constraints <- lapply(problem@constraints, function(c) { canonical_form(c)[[2]] })
  constraints <- unlist(constraints, recursive = TRUE)
  data$objective <- objective
  data$constraints <- constraints
  variables <- variables(problem)[[1]]
  data[BOOL_IDX] <- lapply(variables@boolean_idx, function(t) { t[1] })
  data[INT_IDX] <- lapply(variables@integer_idx, function(t) { t[1] })

  # Order and group constraints.
  inv_data <- list()
  inv_data[var_id(object)] <- id(variables(problem)[[1]])
  inv_data[EQ_CONSTR] <- problem@constraints
  inv_data$is_mip <- is_mixed_integer(problem)
  return(list(data, inv_data))
})

setMethod("invert", signature(object = "XPRESS", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  status <- solution[STATUS]

  primal_vars <- NA
  dual_vars <- NA
  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution[VALUE]
    primal_vars <- list()
    primal_vars[inverse_data[var_id(object)]] <- solution$primal
    if(!is_mip(inverse_data))
      dual_vars <- get_dual_values(solution[EQ_DUAL], extract_dual_value, inverse_data[EQ_CONSTR])
  } else {
    if(status == INFEASIBLE)
      opt_val <- Inf
    else if(status == UNBOUNDED)
      opt_val <- -Inf
    else
      opt_val <- NA
  }

  other <- list()
  other[XPRESS_IIS] <- solution[XPRESS_IIS]
  other[XPRESS_TROW] <- solution[XPRESS_TROW]
  return(Solution(status, opt_val, primal_vars, dual_vars, other))
})

setMethod("solve_via_data", "XPRESS", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  solver <- XPRESS_OLD()
  solver_opts[BOOL_IDX] <- data[BOOL_IDX]
  solver_opts[INT_IDX] <- data[INT_IDX]
  prob_data <- list()
  prob_data[name(object)] <- ProblemData()
  solve(solver, data$objective, data$constraints, prob_data, warm_start, verbose, solver_opts)
})
