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

ConicSolver.get_coeff_offset(expr) {
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

ConicSolver.get_spacing_matrix(shape, spacing, offset) {
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
    eq_dual <- get_dual_values(solution$eq_dual, extract_dual_value, inverse_data[EQ_CONSTR])
    leq_dual <- get_dual_values(solution$ineq_dual, extract_dual_value, inverse_data[NEQ_CONSTR])
    update(eq_dual, leq_dual)
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
  requireNamespace(rcbc, quietly = TRUE)
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
  inv_data[EQ_CONSTR] <- eq_constr
  leq_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "NonPos" })]
  inv_data[NEQ_CONSTR] <- leq_constr
  return(list(data, inv_data))
})

setMethod("invert", signature(object = "CBC", solution = "Solution", inverse_data = "InverseData") {
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
setMethod("import_solver", "CPLEX", function(object) { requireNamespace(Rcplex, quietly = TRUE) })
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
      update(eq_dual, leq_dual)
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
setMethod("import_solver", "CVXOPT", function(object) { requireNamespace(cccopt, quietly = TRUE) })

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
