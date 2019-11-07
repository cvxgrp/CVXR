is_stuffed_cone_constraint <- function(constraint) {
  # Conic solvers require constraints to be stuffed in the following way.
  if(length(variables(constraint)) != 1)
    return(FALSE)
  for(arg in constraint@args) {
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
  expr <- expr(objective)
  return(is_affine(expr) && length(variables(expr)) == 1 && class(expr) == "AddExpression" && length(expr@args) == 2
                         && class(expr@args[[1]]) %in% c("MulExpression", "Multiply") && class(expr@args[[2]]) == "Constant")
}

#'
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
.ConeDims <- setClass("ConeDims", representation(constr_map = "list", zero = "numeric", nonpos = "numeric", exp = "numeric", soc = "list", psd = "list"),
                                  prototype(zero = NA_real_, nonpos = NA_real_, exp = NA_real_, soc = list(), psd = list()))

ConeDims <- function(constr_map) { .ConeDims(constr_map = constr_map) }

setMethod("initialize", "ConeDims", function(.Object, constr_map, zero = NA_real_, nonpos = NA_real_, exp = NA_real_, soc = list(), psd = list()) {
  .Object@zero <- ifelse(is.null(constr_map$ZeroConstraint), 0, sum(sapply(constr_map$ZeroConstraint, function(c) { size(c) })))
  .Object@nonpos <- ifelse(is.null(constr_map$NonPosConstraint), 0, sum(sapply(constr_map$NonPosConstraint, function(c) { size(c) })))
  .Object@exp <- ifelse(is.null(constr_map$ExpCone), 0, sum(sapply(constr_map$ExpCone, function(c) { num_cones(c) })))
  .Object@soc <- as.list(Reduce(c, lapply(constr_map$SOC, function(c) { cone_sizes(c) })))
  .Object@psd <- lapply(constr_map$PSDConstraint, function(c) { dim(c)[1] })
  return(.Object)
})

#'
#' The ConicSolver class.
#'
#' Conic solver class with reduction semantics.
#'
#' @rdname ConicSolver-class
ConicSolver <- setClass("ConicSolver", representation(dims = "character"),   # The key that maps to ConeDims in the data returned by perform().
                                       prototype(dims = "dims"), contains = "ReductionSolver")

# Every conic solver must support Zero and NonPos constraints.
setMethod("supported_constraints", "ConicSolver", function(solver) { c("ZeroConstraint", "NonPosConstraint") })

# Some solvers cannot solve problems that do not have constraints.
# For such solvers, requires_constr should return TRUE.
setMethod("requires_constr", "ConicSolver", function(solver) { FALSE })

setMethod("accepts", signature(object = "ConicSolver", problem = "Problem"), function(object, problem) {
  return(class(problem@objective) == "Minimize" && (mip_capable(object) || !is_mixed_integer(problem)) && is_stuffed_cone_objective(problem@objective)
    && length(convex_attributes(variables(problem))) == 0 && (length(problem@constraints) > 0 || !requires_constr(object))
    && all(sapply(problem@constraints, function(c) { class(c) %in% supported_constraints(object) }))
    && all(sapply(problem@constraints, is_stuffed_cone_constraint)))
})

ConicSolver.get_coeff_offset <- function(expr) {
  # Return the coefficient and offset in A %*% x + b.
  if(class(expr) == "Reshape")   # May be a Reshape as root.
    expr <- expr@args[[1]]
  if(length(expr@args[[1]]@args) == 0) {   # Convert data to float64.
    # expr is t(c) %*% x
    offset <- 0
    coeff <- value(expr@args[[1]])
  } else {
    # expr is t(c) %*% x + d
    offset <- matrix(t(value(expr@args[[2]])), ncol = 1)
    coeff <- value(expr@args[[1]]@args[[1]])
  }
  # Convert scalars to sparse matrices.
  if(is.atomic(coeff) && length(coeff) == 1)
    coeff <- Matrix(coeff, sparse = TRUE)
  return(list(coeff, offset))
}

ConicSolver.get_spacing_matrix <- function(dim, spacing, offset) {
  # Returns a sparse matrix that spaces out an expression.
  val_arr <- c()
  row_arr <- c()
  col_arr <- c()

  # Selects from each column.
  for(var_row in 1:dim[2]) {
    val_arr <- c(val_arr, 1.0)
    row_arr <- c(row_arr, spacing*(var_row - 1) + offset + 1)
    col_arr <- c(col_arr, var_row)
  }

  # Pad out number of rows with zeroes.
  mat <- sparseMatrix(i = row_arr, j = col_arr, x = val_arr)
  if(dim[1] > nrow(mat)) {
    pad <- sparseMatrix(i = c(), j = c(), dims = c(dim[1] - nrow(mat), ncol(mat)))
    mat <- rbind(mat, pad)
  }
  if(!all(dim(mat) == dim))
    stop("Dimension mismatch")
  return(mat)
}

setMethod("reduction_format_constr", "ConicSolver", function(object, problem, constr, exp_cone_order) {
  coeffs <- list()
  offsets <- list()
  for(arg in constr@args) {
    res <- ConicSolver.get_coeff_offset(arg)
    coeffs <- c(coeffs, list(res[[1]]))
    offsets <- c(offsets, list(res[[2]]))
  }
  height <- ifelse(length(coeffs) == 0, 0, sum(sapply(coeffs, function(c) { dim(c)[1] })))

  if(class(constr) %in% c("NonPosConstraint", "ZeroConstraint"))
    # Both of these constraints have but a single argument.
    # t(c) %*% x + b (<)= 0 if and only if t(c) %*% x (<)= b.
    return(list(Matrix(coeffs[[1]], sparse = TRUE), -offsets[[1]]))
  else if(class(constr) == "SOC") {
    # Group each t row with appropriate X rows.
    if(constr@axis != 2)
      stop("SOC must be applied with axis == 2")

    # Interleave the rows of coeffs[[1]] and coeffs[[2]]:
    #   coeffs[[1]][1,]
    #   coeffs[[2]][1:(gap-1),]
    #   coeffs[[1]][2,]
    #   coeffs[[2]][gap:(2*(gap-1)),]
    # TODO: Keep X_coeff sparse while reshaping!
    X_coeff <- coeffs[[2]]
    reshaped <- matrix(t(X_coeff), nrow = nrow(coeffs[[1]]), byrow = TRUE)
    stacked <- -cbind(coeffs[[1]], reshaped)
    stacked <- matrix(t(stacked), nrow = nrow(coeffs[[1]]) + nrow(X_coeff), ncol = ncol(coeffs[[1]]), byrow = TRUE)

    offset <- cbind(offsets[[1]], matrix(t(offsets[[2]]), nrow = nrow(offsets[[1]]), byrow = TRUE))
    offset <- matrix(t(offset), ncol = 1)
    return(list(Matrix(stacked, sparse = TRUE), as.vector(offset)))
  } else if(class(constr) == "ExpCone") {
    for(i in 1:length(coeffs)) {
      mat <- ConicSolver.get_spacing_matrix(c(height, nrow(coeffs[[i]])), length(exp_cone_order), exp_cone_order[i])
      offsets[[i]] <- mat %*% offsets[[i]]
      coeffs[[i]] <- -mat %*% coeffs[[i]]
    }
    # return(list(sum(coeffs), sum(offsets)))
    coeff_sum <- Reduce("+", coeffs)
    offset_sum <- Reduce("+", offsets)
    return(list(coeff_sum, as.vector(offset_sum)))
  } else
    # subclasses must handle PSD constraints.
    stop("Unsupported constraint type.")
})

setMethod("group_coeff_offset", "ConicSolver", function(object, problem, constraints, exp_cone_order) {
  # Combine the constraints into a single matrix, offset.
  if(is.na(constraints) || is.null(constraints) || length(constraints) == 0)
    return(list(Matrix(0, nrow = 0, ncol = 0), numeric(0)))

  matrices <- list()
  offset <- c()
  for(cons in constraints) {
    res <- reduction_format_constr(object, problem, cons, exp_cone_order)
    matrices <- c(matrices, list(res[[1]]))
    offset <- c(offset, res[[2]])
  }
  coeff <- Matrix(do.call(rbind, matrices), sparse = TRUE)
  return(list(coeff, offset))
})

setMethod("invert", signature(object = "ConicSolver", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  # Returns the solution to the original problem given the inverse_data.
  status <- solution$status

  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution$value
    primal_vars <- list()
    primal_vars[[inverse_data[[object@var_id]]]] <- solution$primal
    eq_dual <- get_dual_values(solution$eq_dual, extract_dual_value, inverse_data[object@eq_constr])
    leq_dual <- get_dual_values(solution$ineq_dual, extract_dual_value, inverse_data[object@neq_constr])
    eq_dual <- utils::modifyList(eq_dual, leq_dual)
    dual_vars <- eq_dual
  } else {
    primal_vars <- list()
    primal_vars[[inverse_data[[object@var_id]]]] <- NA_real_
    dual_vars <- NA

    if(status == INFEASIBLE)
      opt_val <- Inf
    else if(status == UNBOUNDED)
      opt_val <- -Inf
    else
      opt_val <- NA_real_
  }
  return(Solution(status, opt_val, primal_vars, dual_vars, list()))
})

ECOS <- setClass("ECOS", representation(exp_cone_order = "numeric"),   # Order of exponential cone arguments for solver. Internal only!
                 prototype(exp_cone_order = c(0, 2, 1)), contains = "ConicSolver")

# Solver capabilities.
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

setMethod("import_solver", "ECOS", function(solver) { requireNamespace("ECOSolveR", quietly = TRUE) })

setMethod("name", "ECOS", function(x) { ECOS_NAME })
setMethod("perform", signature(object = "ECOS", problem = "Problem"), function(object, problem) {
  data <- list()
  inv_data <- list()
  inv_data[[object@var_id]] <- id(variables(problem)[[1]])
  offsets <- ConicSolver.get_coeff_offset(problem@objective@args[[1]])
  data[[C_KEY]] <- as.vector(offsets[[1]])
  data[[OFFSET]] <- offsets[[2]]
  inv_data[[OFFSET]] <- data[[OFFSET]][1]

  constr_map <- group_constraints(problem@constraints)
  data[[ConicSolver()@dims]] <- ConeDims(constr_map)

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

setMethod("solve_via_data", "ECOS", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
  requireNamespace("ECOSolveR", quietly = TRUE)
  cones <- ECOS.dims_to_solver_dict(data[[ConicSolver()@dims]])
  ecos_opts <- ECOSolveR::ecos.control()
  ecos_opts$VERBOSE <- as.integer(verbose)
  ecos_opts[names(solver_opts)] <- solver_opts
  solution <- ECOSolveR::ECOS_csolve(c = data[[C_KEY]], G = data[[G_KEY]], h = data[[H_KEY]], dims = cones, A = data[[A_KEY]], b = data[[B_KEY]], control = ecos_opts)
  return(solution)
})

SCS <- setClass("SCS", representation(exp_cone_order = "numeric"),   # Order of exponential cone arguments for solver. Internal only!
                prototype(exp_cone_order = c(0, 1, 2)), contains = "ConicSolver")

# Solver capabilities.
setMethod("mip_capable", "SCS", function(solver) { FALSE })
setMethod("requires_constr", "SCS", function(solver) { TRUE })
setMethod("supported_constraints", "SCS", function(solver) { c(supported_constraints(ConicSolver()), "SOC", "ExpCone", "PSDConstraint") })

# Map of SCS status to CVXR status.
setMethod("status_map", "SCS", function(solver, status) {
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

setMethod("name", "SCS", function(x) { SCS_NAME })
setMethod("import_solver", "SCS", function(solver) { requireNamespace("scs", quietly = TRUE) })

setMethod("reduction_format_constr", "SCS", function(object, problem, constr, exp_cone_order) {
  # Extract coefficient and offset vector from constraint.
  # Special cases PSD constraints, as SCS expects constraints to be
  # imposed on solely the lower triangular part of the variable matrix.
  # Moreover, it requires the off-diagonal coefficients to be scaled by
  # sqrt(2).
  if(is(constr, "PSDConstraint")) {
    expr <- expr(constr)
    triangularized_expr <- scaled_lower_tri(expr + t(expr))/2
    extractor <- CoeffExtractor(InverseData(problem))
    Ab <- affine(extractor, triangularized_expr)
    A_prime <- Ab[[1]]
    b_prime <- Ab[[2]]

    # SCS requests constraints to be formatted as Ax + s = b,
    # where s is constrained to reside in some cone. Here, however,
    # we are formatting the constraint as A"x + b" = -Ax + b; h ence,
    # A = -A", b = b".
    return(list(-1*A_prime, b_prime))
  } else
    callNextMethod(object, problem, constr, exp_cone_order)
})

setMethod("perform", signature(object = "SCS", problem = "Problem"), function(object, problem) {
  # Returns a new problem and data for inverting the new solution.
  data <- list()
  inv_data <- list()
  inv_data[[object@var_id]] <- id(variables(problem)[[1]])

  # Parse the coefficient vector from the objective.
  offsets <- ConicSolver.get_coeff_offset(problem@objective@args[[1]])
  data[[C_KEY]] <- offsets[[1]]
  data[[OFFSET]] <- offsets[[2]]
  data[[C_KEY]] <- as.vector(data[[C_KEY]])
  inv_data[[OFFSET]] <- data[[OFFSET]][1]

  # Order and group nonlinear constraints.
  constr_map <- group_constraints(problem@constraints)
  data[[ConicSolver()@dims]] <- ConeDims(constr_map)
  inv_data[[ConicSolver()@dims]] <- data[[ConicSolver()@dims]]

  # SCS requires constraints to be specified in the following order:
  # 1) Zero cone.
  # 2) Non-negative orthant.
  # 3) SOC.
  # 4) PSD.
  # 5) Exponential.
  zero_constr <- constr_map$ZeroConstraint
  neq_constr <- c(constr_map$NonPosConstraint, constr_map$SOC, constr_map$PSDConstraint, constr_map$ExpCone)
  inv_data[[object@eq_constr]] <- zero_constr
  inv_data[[object@neq_constr]] <- neq_constr

  # Obtain A, b such that Ax + s = b, s \in cones.
  # Note that SCS mandates that the cones MUST be ordered with
  # zero cones first, then non-negative orthant, then SOC, then
  # PSD, then exponential.
  offsets <- group_coeff_offset(object, problem, c(zero_constr, neq_constr), object@exp_cone_order)
  data[[A_KEY]] <- offsets[[1]]
  data[[B_KEY]] <- offsets[[2]]
  return(list(object, data, inv_data))
})

SCS.extract_dual_value <- function(result_vec, offset, constraint) {
  # Extracts the dual value for constraint starting at offset.
  # Special cases PSD constraints, as per the SCS specification.
  if(is(constraint, "PSDConstraint")) {
    dim <- nrow(constraint)
    lower_tri_dim <- floor(dim*(dim+1)/2)
    new_offset <- offset + lower_tri_dim
    lower_tri <- result_vec[(offset + 1):new_offset]
    full <- tri_to_full(lower_tri, dim)
    return(list(full, new_offset))
  } else
    return(extract_dual_value(result_vec, offset, constraint))
}

setMethod("invert", signature(object = "SCS", solution = "list", inverse_data = "list"), function(object, solution, inverse_data) {
  # Returns the solution to the original problem given the inverse_data.
  status <- status_map(object, solution$info$status)

  attr <- list()
  attr[[SOLVE_TIME]] <- solution$info$solveTime
  attr[[SETUP_TIME]] <- solution$info$setupTime
  attr[[NUM_ITERS]] <- solution$info$iter

  if(status %in% SOLUTION_PRESENT) {
    primal_val <- solution$info$pobj
    opt_val <- primal_val + inverse_data[[OFFSET]]
    primal_vars <- list()
    var_id <- inverse_data[[object@var_id]]
    primal_vars[[as.character(var_id)]] <- as.matrix(solution$x)

    num_zero <- inverse_data[[ConicSolver()@dims]]@zero
    eq_idx <- seq_len(num_zero)
    ineq_idx <- seq(num_zero + 1, length.out = length(solution$y) - num_zero)
    eq_dual_vars <- get_dual_values(solution$y[eq_idx], SCS.extract_dual_value, inverse_data[[object@eq_constr]])
    ineq_dual_vars <- get_dual_values(solution$y[ineq_idx], SCS.extract_dual_value, inverse_data[[object@neq_constr]])

    dual_vars <- list()
    dual_vars <- utils::modifyList(dual_vars, eq_dual_vars)
    dual_vars <- utils::modifyList(dual_vars, ineq_dual_vars)
    return(Solution(status, opt_val, primal_vars, dual_vars, attr))
  } else
    return(failure_solution(status))
})

setMethod("solve_via_data", "SCS", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
  requireNamespace("scs", quietly = TRUE)

  # TODO: Cast A to dense because scs R package rejects sparse matrices?
  args <- list(A = data[[A_KEY]], b = data[[B_KEY]], c = data[[C_KEY]])
  if(warm_start && !is.null(solver_cache) && length(solver_cache) > 0 && name(object) %in% names(solver_cache)) {
    args$x <- solver_cache[[name(object)]]$x
    args$y <- solver_cache[[name(object)]]$y
    args$s <- solver_cache[[name(object)]]$s
  }
  cones <- SCS.dims_to_solver_dict(data[[ConicSolver()@dims]])

  # Default to eps = 1e-4 instead of 1e-3.
  if(is.null(solver_opts$eps))
    solver_opts$eps <- 1e-4
  if(!missing(verbose))
    solver_opts$verbose <- verbose

  # Default to acceleration_lookback = 10L instead of (Naras' mistake) 20L!
  # REMOVE when scs R package goes back to correct default
  if(is.null(solver_opts$acceleration_lookback))
    solver_opts$acceleration_lookback <- 10L

  if(!missing(verbose))
    solver_opts$verbose <- verbose

  # Returns the result of the call to the solver.
  results <- scs::scs(A = args$A, b = args$b, obj = args$c, cone = cones, control = solver_opts)
  if(!is.null(solver_cache) && length(solver_cache) > 0)
    solver_cache[[name(object)]] <- results
  return(results)
})

CBC_CONIC <- setClass("CBC_CONIC", contains = "SCS")

# Solver capabilities.
setMethod("mip_capable", "CBC_CONIC", function(solver) { TRUE })

#DK: CBC_Conic to CVXR status. Check these are correct
setMethod("status_map", "CBC_CONIC", function(solver, status) {
  if(status$is_proven_optimal)
    OPTIMAL
  else if(status$is_proven_dual_infeasible || status$is_proven_infeasible)
    INFEASIBLE
  else
    SOLVER_ERROR # probably need to check this the most
})

# Map of CBC_CONIC MIP/LP status to CVXR status.
setMethod("status_map_mip", "CBC_CONIC", function(solver, status) {
  if(status == "solution")
    OPTIMAL
  else if(status == "relaxation_infeasible")
    INFEASIBLE
  else if(status == "stopped_on_user_event")
    SOLVER_ERROR
  else
    stop("CBC_CONIC MIP status unrecognized: ", status)
})

setMethod("status_map_lp", "CBC_CONIC", function(solver, status) {
  if(status == "optimal")
    OPTIMAL
  else if(status == "primal_infeasible")
    INFEASIBLE
  else if(status == "stopped_due_to_errors" || status == "stopped_by_event_handler")
    SOLVER_ERROR
  else
    stop("CBC_CONIC LP status unrecognized: ", status)
})

setMethod("name", "CBC_CONIC", function(x) { CBC_NAME })
setMethod("import_solver", "CBC_CONIC", function(solver) {
  installed <- requireNamespace("rcbc", quietly = TRUE)
  if(!installed)
    stop("Required R package rcbc not found. Please install from https://github.com/dirkschumacher/rcbc")
  return(installed)
})

setMethod("accepts", signature(object = "CBC_CONIC", problem = "Problem"), function(object, problem) {
  # Can CBC_CONIC solve the problem?
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

setMethod("perform", signature(object = "CBC_CONIC", problem = "Problem"), function(object, problem) {
  tmp <- callNextMethod(object, problem)
  object <- tmp[[1]]
  data <- tmp[[2]]
  inv_data <- tmp[[3]]
  variables <- variables(problem)[[1]]
  data[[BOOL_IDX]] <- lapply(variables@boolean_idx, function(t) { t[1] })
  data[[INT_IDX]] <- lapply(variables@integer_idx, function(t) { t[1] })
  inv_data$is_mip <- length(data[[BOOL_IDX]]) > 0 || length(data[[INT_IDX]]) > 0
  return(list(object, data, inv_data))
})

#DK: Changed solution class from "Solution" to rcbc_milp_result. Same with CBC_Conic  to
#setMethod("invert", "CBC_CONIC",  function(object, solution, inverse_data) {
setMethod("invert", signature(object = "CBC_CONIC", solution = "list", inverse_data = "list"),  function(object, solution, inverse_data) {
  # Returns the solution to the original problem given the inverse_data.
  solution <- solution[[1]]
  status <- status_map(object, solution)

  primal_vars <- list()
  dual_vars <- list()
  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution$objective_value
    primal_vars[[object@var_id]] <- solution$column_solution
  } else {
    if(status == INFEASIBLE)
      opt_val <- Inf
    else if(status == UNBOUNDED)
      opt_val <- -Inf
    else
      opt_val <- NA_real_
  }

  return(Solution(status, opt_val, primal_vars, dual_vars, list()))
})

setMethod("solve_via_data", "CBC_CONIC", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
  requireNamespace("rcbc", quietly = TRUE)

  cvar <- data$c
  b <- data$b
  A <- data$A
  dims <- SCS.dims_to_solver_dict(data$dims)

  if(is.null(dim(data$c))){
    n <- length(cvar) # Should dim be used here?
  } else {
    n <- dim(cvar)[1]
  }

  # Initialize variable constraints
  var_lb <- rep(-Inf, n)
  var_ub <- rep(Inf, n)
  is_integer <- rep.int(FALSE, n)
  row_ub <- rep(Inf, nrow(A))
  row_lb <- rep(-Inf, nrow(A))

  #Setting equality constraints
  if(dims[[EQ_DIM]] > 0){
    row_ub[1:dims[[EQ_DIM]]] <- b[1:dims[[EQ_DIM]]]
    row_lb[1:dims[[EQ_DIM]]] <- b[1:dims[[EQ_DIM]]]
  }

  #Setting inequality constraints
  leq_start <- dims[[EQ_DIM]]
  leq_end <- dims[[EQ_DIM]] + dims[[LEQ_DIM]]
  if(leq_start != leq_end){
    row_ub[(leq_start+1):(leq_end)] <- b[(leq_start+1):(leq_end)]
  }

  # Make boolean constraints
  if(length(data$bool_vars_idx) > 0){
    var_lb[unlist(data$bool_vars_idx)] <- 0
    var_ub[unlist(data$bool_vars_idx)] <- 1
    is_integer[unlist(data$bool_vars_idx)] <- TRUE
  }

  if(length(data$int_vars_idx) > 0) {
    is_integer[unlist(data$int_vars_idx)] <- TRUE
  }

  result <- rcbc::cbc_solve(
    obj = cvar,
    mat = A,
    row_ub = row_ub,
    row_lb = row_lb,
    col_lb = var_lb,
    col_ub = var_ub,
    is_integer = is_integer,
    max = FALSE
  )

  # solver <- 'blah'
  # solver_opts[[BOOL_IDX]] <- data[[BOOL_IDX]]
  # solver_opts[[INT_IDX]] <- data[[INT_IDX]]
  # prob_data <- list()
  # prob_data[[name(object)]] <- ProblemData()
  # return(solve(solver, data$objective, data$constraints, prob_data, warm_start, verbose, solver_opts))

  return(list(result))
})

CPLEX_CONIC <- setClass("CPLEX_CONIC", contains = "SCS")

setMethod("mip_capable", "CPLEX_CONIC", function(solver) { TRUE })
setMethod("supported_constraints", "CPLEX_CONIC", function(solver) { c(supported_constraints(ConicSolver()), "SOC") })
setMethod("name", "CPLEX_CONIC", function(x) { CPLEX_NAME })
setMethod("import_solver", "CPLEX_CONIC", function(solver) { requireNamespace("Rcplex", quietly = TRUE) })
setMethod("accepts", signature(object = "CPLEX_CONIC", problem = "Problem"), function(object, problem) {
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

# Map of CPLEX status to CVXR status.
# TODO: Add more!
setMethod("status_map", "CPLEX_CONIC", function(solver, status) {
  if(status %in% c(1, 101)){
    OPTIMAL
  } else if(status %in% c(3, 22, 4, 103)){
    INFEASIBLE
  } else if(status %in% c(2, 21, 118)){
    UNBOUNDED
  } else if(status %in% c(10, 107)){
    USER_LIMIT
  } else
    stop("CPLEX status unrecognized: ", status)
})

setMethod("perform", signature(object = "CPLEX_CONIC", problem = "Problem"), function(object, problem) {
  #COPIED OVER FROM SCS CONIC, which is what CVXPY does (except they superclass it to a class w the same method)

  # Returns a new problem and data for inverting the new solution.
  data <- list()
  inv_data <- list()
  inv_data[[object@var_id]] <- id(variables(problem)[[1]])

  # Parse the coefficient vector from the objective.
  offsets <- ConicSolver.get_coeff_offset(problem@objective@args[[1]])
  data[[C_KEY]] <- offsets[[1]]
  data[[OFFSET]] <- offsets[[2]]
  data[[C_KEY]] <- as.vector(data[[C_KEY]])
  inv_data[[OFFSET]] <- data[[OFFSET]][1]

  # Order and group nonlinear constraints.
  constr_map <- group_constraints(problem@constraints)
  data[[ConicSolver()@dims]] <- ConeDims(constr_map)
  inv_data[[ConicSolver()@dims]] <- data[[ConicSolver()@dims]]

  # SCS requires constraints to be specified in the following order:
  # 1) Zero cone.
  # 2) Non-negative orthant.
  # 3) SOC.
  # 4) PSD.
  # 5) Exponential.
  zero_constr <- constr_map$ZeroConstraint
  neq_constr <- c(constr_map$NonPosConstraint, constr_map$SOC, constr_map$PSDConstraint, constr_map$ExpCone)
  inv_data[[object@eq_constr]] <- zero_constr
  inv_data[[object@neq_constr]] <- neq_constr
  inv_data$is_mip <- length(data[[BOOL_IDX]]) > 0 || length(data[[INT_IDX]]) > 0

  # Obtain A, b such that Ax + s = b, s \in cones.
  # Note that SCS mandates that the cones MUST be ordered with
  # zero cones first, then non-negative orthant, then SOC, then
  # PSD, then exponential.
  offsets <- group_coeff_offset(object, problem, c(zero_constr, neq_constr), object@exp_cone_order)
  data[[A_KEY]] <- offsets[[1]]
  data[[B_KEY]] <- offsets[[2]]
  return(list(object, data, inv_data))
})

setMethod("invert", signature(object = "CPLEX_CONIC", solution = "list", inverse_data = "list"), function(object, solution, inverse_data) {
  # Returns the solution to the original problem given the inverse_data.
  model <- solution$model

  status <- status_map(object, model$status)

  if(status %in% SOLUTION_PRESENT) {
    #Get objective value
    opt_val <- model$obj + inverse_data[[OFFSET]]

    #Get solution
    primal_vars <- list()
    primal_vars[[as.character(inverse_data[[object@var_id]])]] <- model$xopt

    #Only add duals if not a MIP
    if(!inverse_data[["is_mip"]]) {
      #eq and leq constraints all returned at once by CPLEX
      #eq_dual <- get_dual_values(solution$eq_dual, extract_dual_value, inverse_data[[object@eq_constr]])
      #leq_dual <- get_dual_values(solution$ineq_dual, extract_dual_value, inverse_data[[object@neq_constr]])
      #eq_dual <- utils::modifyList(eq_dual, leq_dual)
      dual_vars <- get_dual_values(solution$y, extract_dual_value, inverse_data[[object@neq_constr]])
      #dual_vars <- eq_dual

    } else {
      primal_vars <- list()
      primal_vars[[inverse_data[[object@var_id]]]] <- NA_real_
      if(!inverse_data[["is_mip"]]) {
        dual_var_ids <- sapply(c(inverse_data[[object@eq_constr]], inverse_data[[object@neq_constr]]), function(constr) { constr@id })
        dual_vars <- as.list(rep(NA_real_, length(dual_var_ids)))
        names(dual_vars) <- dual_var_ids
      }

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

setMethod("solve_via_data", "CPLEX_CONIC", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
  cvar <- data[[C_KEY]]
  bvec <- data[[B_KEY]]
  Amat <- data[[A_KEY]]
  dims <- data[[DIMS]]

  total_soc <- sum(unlist(dims@soc))
  n_var <- length(cvar)
  cvar <- c(cvar, rep(0, total_soc))

  #Initializing variable types
  vtype <- rep("C", n_var + total_soc)

  #Setting Boolean variable types
  for(i in seq_along(data[BOOL_IDX]$bool_vars_idx)){
    vtype[data[BOOL_IDX]$bool_vars_idx[i]] <- "B"
  }
  #Setting Integer variable types
  for(i in seq_along(data[INT_IDX]$int_vars_idx)){
    vtype[data[BOOL_IDX]$int_vars_idx[i]] <- "I"
  }

  #Setting sense of the A matrix
  sense_vec <- rep("E", nrow(Amat))

  #Add inequalities
  leq_start <- dims@zero

  for(j in 1:dims@nonpos){
    sense_vec[leq_start + j] <- "L"
  }

  #Setting Lower bounds of variables
  lb <- rep(-Inf, n_var + total_soc)

  qc <- list()

  soc_start <- leq_start + dims@nonpos
  current_var <- n_var

  for(i in seq_along(dims@soc)){
    for(j in 1:dims@soc[[i]]){
      sense_vec[soc_start + j] <- "E"
      if(j == 1){
        lb[current_var + j] <- 0 #The first variable of every SOC has a 0 lower bound
      }
    }

    #Add SOC vars to linear constraints
    n_soc <- dims@soc[[i]]
    Asub <- matrix(0, nrow = nrow(Amat), ncol = n_soc)
    Asub[(soc_start+1):(soc_start + n_soc),] <- diag(rep(1, n_soc))
    Amat <- cbind(Amat, Asub)

    #Add quadratic constraints
    qc_mat <- matrix(0, nrow = n_var + total_soc, ncol = n_var + total_soc)
    qc_mat[current_var + 1, current_var + 1] <- -1
    for(k in 1:(n_soc-1)){
      qc_mat[current_var + 1 + k, current_var + 1 + k] <- 1
    }
    qc[[i]] <- qc_mat

    soc_start <- soc_start + n_soc
    current_var <- current_var + n_soc
  }

  QC <- list(QC = list(Q = qc), dir = rep("L", length(dims@soc)) , b = rep(0.0, length(dims@soc)))

  #Setting verbosity off
  if(!verbose){
    control <- list(trace=0)
  } else{
    control <- list(trace=1)
  }

  #David: not sure what to do with this part
  # TODO: The code in CVXR/problems/solvers.R sets CPLEX parameters in the same way,
  # and the code is duplicated here. This should be refactored.
  if(length(solver_opts) > 0){
    kwargs <- sort(names(solver_opts))
    if("cplex_params" %in% kwargs) {
      for(param in names(solver_opts$cplex_params)) {
        value <- solver_opts$cplex_params[param]
        tryCatch({
          eval(paste("set(model@parameters@", param, ", value)", sep = ""))
        }, error = function(e) {
          stop("Invalid CPLEX parameter, value pair (", param, ", ", value, ")")
        })
      }
      kwargs$cplex_params <- NULL
    }

    if("cplex_filename" %in% kwargs) {
      filename <- solver_opts$cplex_filename
      if(!is.na(filename) && !is.null(filename))
        write(model, filename)
      kwargs$cplex_filename <- NULL
    }

    if(length(is.null(kwargs))<=0 || length(is.na(kwargs))<=0)
      stop("Invalid keyword argument ", kwargs[[1]])
  }

  # Solve problem.
  results_dict <- list()

  tryCatch({
    # Define CPLEX problem and solve
    model <- Rcplex::Rcplex_solve_QCP(cvec=cvar, Amat=Amat, bvec=bvec, QC=QC, lb=lb, ub=Inf,
                            sense=sense_vec, objsense="min", vtype=vtype, control=control)
    #control parameter would be used to set specific solver arguments. See cran Rcplex documentation
  }, error = function(e) {
    results_dict$status <- SOLVER_ERROR
  }
  )

  #Changing dualvar to include SOC
  y <- model$extra$lambda
  soc_start <- leq_start + dims@nonpos
  for(i in seq_along(dims@soc)){
    y <- append(y, 0, soc_start)
    soc_start <- soc_start + dims@soc[[i]] + 1
  }
  results_dict$y <- -y
  #results_dict$eq_dual <- y[1:dims@zero]
  #results_dict$ineq_dual <- y[-(1:dims@zero)]
  results_dict$model <- model

  return(results_dict)

})

setClass("CVXOPT", contains = "ECOS")
CVXOPT <- function() { new("CVXOPT") }

# Solver capabilities.
setMethod("mip_capable", "CVXOPT", function(solver) { FALSE })
setMethod("supported_constraints", "CVXOPT", function(solver) { c(supported_constraints(ConicSolver()), "SOC", "ExpCone", "PSDConstraint") })

# Map of CVXOPT status to CVXR status.
setMethod("status_map", "CVXOPT", function(solver, status) {
  if(status == "optimal")
    OPTIMAL
  else if(status == "infeasible")
    INFEASIBLE
  else if(status == "unbounded")
    UNBOUNDED
  else if(status == "solver_error")
    SOLVER_ERROR
  else
    stop("CVXOPT status unrecognized: ", status)
})

setMethod("name", "CVXOPT", function(x) { CVXOPT_NAME })
setMethod("import_solver", "CVXOPT", function(solver) { requireNamespace("cccopt", quietly = TRUE) })

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

setMethod("perform", signature(object = "CVXOPT", problem = "Problem"), function(object, problem) {
  data <- list()
  inv_data <- list()
  inv_data[[object@var_id]] <- id(variables(problem)[[1]])
  tmp <- ConicSolver.get_coeff_offset(problem@objective@args[[1]])
  data[[C_KEY]] <- as.vector(tmp[[1]])
  data[[OFFSET]] <- tmp[[2]]
  inv_data[[OFFSET]] <- data[[OFFSET]][1]

  constr_map <- group_constraints(problem@constraints)
  data[[ConicSolver()@dims]] <- ConeDims(constr_map)

  inv_data[[object@eq_constr]] <- constr_map$ZeroConstraint
  tmp <- group_coeff_offset(object, problem, constr_map$ZeroConstraint, ECOS()@exp_cone_order)
  data[[A_KEY]] <- tmp[[1]]
  data[[B_KEY]] <- tmp[[2]]

  # Order and group nonlinear constraints.
  neq_constr <- c(constr_map$NonPosConstraint, constr_map$SOC, constr_map$PSDConstraint)
  inv_data[[object@neq_constr]] <- neq_constr
  tmp <- group_coeff_offset(object, problem, neq_constr, ECOS()@exp_cone_order)
  data[[G_KEY]] <- tmp[[1]]
  data[[H_KEY]] <- tmp[[2]]

  var <- variables(problem)[[1]]
  data[[BOOL_IDX]] <- as.integer(var@boolean_idx[,1])
  data[[INT_IDX]] <- as.integer(var@integer_idx[,1])

  #Add information about
  return(list(object, data, inv_data))
})

setMethod("solve_via_data", "CVXOPT", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
  solver <- CVXOPT_OLD()
  prob_data <- list()
  prob_data[[name(object)]] <- ProblemData()
  solve(solver, data$objective, data$constraints, prob_data, warm_start, verbose, solver_opts)
})

ECOS_BB <- setClass("ECOS_BB", contains = "ECOS")

setMethod("mip_capable", "ECOS_BB", function(solver) { TRUE })
setMethod("name", "ECOS_BB", function(x) { ECOS_BB_NAME })
setMethod("perform", signature(object = "ECOS_BB", problem = "Problem"), function(object, problem) {
  res <- callNextMethod(object, problem)
  object <- res[[1]]
  data <- res[[2]]
  inv_data <- res[[3]]

  # Because the problem variable is single dimensional, every
  # boolean/integer index has length one.
  var <- variables(problem)[[1]]
  data[[BOOL_IDX]] <- as.integer(var@boolean_idx[,1])
  data[[INT_IDX]] <- as.integer(var@integer_idx[,1])
  return(list(object, data, inv_data))
})

setMethod("solve_via_data", "ECOS_BB", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
  requireNamespace("ECOSolveR", quietly = TRUE)
  cones <- ECOS.dims_to_solver_dict(data[[ConicSolver()@dims]])
  ecos_opts <- ECOSolveR::ecos.control()
  ecos_opts$VERBOSE <- as.integer(verbose)
  ecos_opts[names(solver_opts)] <- solver_opts
  solution <- ECOSolveR::ECOS_csolve(c = data[[C_KEY]], G = data[[G_KEY]], h = data[[H_KEY]], dims = cones, A = data[[A_KEY]], b = data[[B_KEY]],
                                     bool_vars = data[[BOOL_IDX]], int_vars = data[[INT_IDX]], control = ecos_opts)
  return(solution)
})

# Utility method for formatting a ConeDims instance into a dictionary
# that can be supplied to ECOS.
ECOS.dims_to_solver_dict <- function(cone_dims) {
  cones <- list(l = as.integer(cone_dims@nonpos),
                q = lapply(cone_dims@soc, function(v) { as.integer(v) }),
                e = as.integer(cone_dims@exp))
  return(cones)
}

GLPK <- setClass("GLPK", contains = "CVXOPT")
setMethod("mip_capable", "GLPK", function(solver) { FALSE })
setMethod("supported_constraints", "GLPK", function(solver) { supported_constraints(ConicSolver()) })

# Map of GLPK_MI status to CVXR status.
setMethod("status_map", "GLPK", function(solver, status) {
  if(status == 5)
    OPTIMAL
  else if(status == 2)
    SOLUTION_PRESENT
  else if(status == 3 | status == 4)
    INFEASIBLE
  else if(status == 1 | status == 6)
    UNBOUNDED
  else
    stop("GLPK status unrecognized: ", status)
})

setMethod("name", "GLPK", function(x) { GLPK_NAME })
setMethod("import_solver", "GLPK", function(solver) { requireNamespace("Rglpk", quietly = TRUE) })

setMethod("invert", signature(object = "GLPK", solution = "list", inverse_data = "list"), function(object, solution, inverse_data) {
  status <- solution$status
  primal_vars <- list()
  dual_vars <- list()
  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution$value
    primal_vars <- list()
    primal_vars[[as.character(inverse_data[[as.character(object@var_id)]])]] <- solution$primal
    return(Solution(status, opt_val, primal_vars, dual_vars, list()))
  } else
    return(failure_solution(status))
})

setMethod("solve_via_data", "GLPK", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
  if(verbose)
    solver_opts$verbose <- verbose
  solver_opts$canonicalize_status <- FALSE

  # Construct problem data.
  c <- data[[C_KEY]]
  dims <- data[[ConicSolver()@dims]]
  nvar <- length(c)
  A <- data[[A_KEY]]
  b <- data[[B_KEY]]
  if(nrow(A) == 0)
    A <- Matrix(0, nrow = 0, ncol = length(c))

  G <- data[[G_KEY]]
  h <- data[[H_KEY]]
  if(nrow(G) == 0)
    G <- Matrix(0, nrow = 0, ncol = length(c))

  mat <- rbind(A, G)
  rhs <- c(b, h)

  bounds <- list(lower = list(ind = seq_along(c), val = rep(-Inf, nvar)))
  types <- rep("C", nvar)
  bools <- data[[BOOL_IDX]]
  ints <- data[[INT_IDX]]


  if (length(bools) > 0) {
    types[bools] <- "B"
  }
  if (length(ints) > 0) {
    types[ints] <- "I"
  }

  results_dict <- Rglpk::Rglpk_solve_LP(obj = c,
                                        mat = slam::as.simple_triplet_matrix(mat),
                                        dir = c(rep("==", dims@zero),
                                                rep("<=", dims@nonpos)),
                                        rhs = rhs,
                                        bounds = bounds,
                                        types = types,
                                        control = solver_opts,
                                        max = FALSE)

  # Convert results to solution format.
  solution <- list()
  solution[[STATUS]] <- status_map(object, results_dict$status)
  if(solution[[STATUS]] %in% SOLUTION_PRESENT) {
    ## Get primal variable values
    solution[[PRIMAL]] <- results_dict$solution
    ## Get objective value
    solution[[VALUE]] <- results_dict$optimum
    # solution[[EQ_DUAL]] <- results_dict$auxiliary[[1]]   # TODO: How do we get the dual variables?
    # solution[[INEQ_DUAL]] <- results_dict$auxiliar[[2]]
  }
  return(solution)
})

GLPK_MI <- setClass("GLPK_MI", contains = "GLPK")
setMethod("mip_capable", "GLPK_MI", function(solver) { TRUE })
setMethod("supported_constraints", "GLPK_MI", function(solver) { supported_constraints(ConicSolver()) })

# Map of GLPK_MI status to CVXR status.
setMethod("status_map", "GLPK_MI", function(solver, status) {
  if(status == 5)
    OPTIMAL
  else if(status == 2)
    SOLUTION_PRESENT #Not sure if feasible is the same thing as this
  else if(status == 3 | status == 4)
    INFEASIBLE
  else if(status == 1 | status == 6)
    UNBOUNDED
  else
    stop("GLPK_MI status unrecognized: ", status)
})

setMethod("name", "GLPK_MI", function(x) { GLPK_MI_NAME })
setMethod("solve_via_data", "GLPK_MI", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
  if(verbose)
    solver_opts$verbose <- verbose
  solver_opts$canonicalize_status <- FALSE

  # Construct problem data.
  c <- data[[C_KEY]]
  dims <- data[[ConicSolver()@dims]]
  nvar <- length(c)
  A <- data[[A_KEY]]
  b <- data[[B_KEY]]
  if(nrow(A) == 0)
    A <- Matrix(0, nrow = 0, ncol = length(c))

  G <- data[[G_KEY]]
  h <- data[[H_KEY]]
  if(nrow(G) == 0)
    G <- Matrix(0, nrow = 0, ncol = length(c))

  mat <- rbind(A, G)
  rhs <- c(b, h)
  bounds <- list(lower = list(ind = seq_along(c), val = rep(-Inf, nvar)))
  types <- rep("C", nvar)
  bools <- data[[BOOL_IDX]]
  ints <- data[[INT_IDX]]
  if (length(bools) > 0) {
    types[bools] <- "B"
  }
  if (length(ints) > 0) {
    types[ints] <- "I"
  }

  results_dict <- Rglpk::Rglpk_solve_LP(obj = c,
                                        mat = slam::as.simple_triplet_matrix(mat),
                                        dir = c(rep("==", dims@zero),
                                                rep("<=", dims@nonpos)),
                                        rhs = rhs,
                                        bounds = bounds,
                                        types = types,
                                        control = solver_opts,
                                        max = FALSE)

  # Convert results to solution format.
  solution <- list()

  solution[[STATUS]] <- status_map(object, results_dict$status)
  if(solution[[STATUS]] %in% SOLUTION_PRESENT) {
    ## Get primal variable values
    solution[[PRIMAL]] <- results_dict$solution
    ## Get objective value
    solution[[VALUE]] <- results_dict$optimum
    # solution[[EQ_DUAL]] <- results_dict$auxiliary[[1]]   # TODO: How do we get the dual variables?
    # solution[[INEQ_DUAL]] <- results_dict$auxiliar[[2]]
    solution[[EQ_DUAL]] <- list()
    solution[[INEQ_DUAL]] <- list()
  }
  solution
})

GUROBI_CONIC <- setClass("GUROBI_CONIC", contains = "SCS")

# Solver capabilities.
setMethod("mip_capable", "GUROBI_CONIC", function(solver) { TRUE })
setMethod("supported_constraints", "GUROBI_CONIC", function(solver) { c(supported_constraints(ConicSolver()), "SOC") })

# Map of Gurobi status to CVXR status.
setMethod("status_map", "GUROBI_CONIC", function(solver, status) {
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

setMethod("name", "GUROBI_CONIC", function(x) { GUROBI_NAME })
setMethod("import_solver", "GUROBI_CONIC", function(solver) { requireNamespace("gurobi", quietly = TRUE) })

# Map of GUROBI status to CVXR status.
setMethod("status_map", "GUROBI_CONIC", function(solver, status) {
  if(status == 2 || status == "OPTIMAL")
    OPTIMAL
  else if(status == 3 || status == 6 || status == "INFEASIBLE") #DK: I added the words because the GUROBI solver seems to return the words
    INFEASIBLE
  else if(status == 5 || status == "UNBOUNDED")
    UNBOUNDED
  else if(status == 4 | status == "INF_OR_UNBD")
    INFEASIBLE_INACCURATE
  else if(status %in% c(7,8,9,10,11,12))
    SOLVER_ERROR   # TODO: Could be anything
  else if(status == 13)
    OPTIMAL_INACCURATE   # Means time expired.
  else
    stop("GUROBI status unrecognized: ", status)
})


setMethod("accepts", signature(object = "GUROBI_CONIC", problem = "Problem"), function(object, problem) {
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

setMethod("perform", signature(object = "GUROBI_CONIC", problem = "Problem"), function(object, problem) {
  tmp <- callNextMethod(object, problem)
  object <- tmp[[1]]
  data <- tmp[[2]]
  inv_data <- tmp[[3]]
  variables <- variables(problem)[[1]]
  data[[BOOL_IDX]] <- lapply(variables@boolean_idx, function(t) { t[1] })
  data[[INT_IDX]] <- lapply(variables@integer_idx, function(t) { t[1] })
  inv_data$is_mip <- length(data[[BOOL_IDX]]) > 0 || length(data[[INT_IDX]]) > 0
  return(list(object, data, inv_data))
})

setMethod("invert", signature(object = "GUROBI_CONIC", solution = "list", inverse_data = "list"), function(object, solution, inverse_data) {

  status <- solution$status
  dual_vars <- list()

  #CVXPY doesn't include for some reason?
  #attr <- list()
  #attr[[SOLVE_TIME]] <- solution$runtime
  #attr[[NUM_ITERS]] <- solution$baritercount

  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution$value + inverse_data[[OFFSET]]
    primal_vars <- list()
    primal_vars[[as.character(inverse_data[[as.character(object@var_id)]])]] <- solution$primal
    if(!inverse_data[["is_mip"]]) {
      eq_dual <- get_dual_values(solution$eq_dual, extract_dual_value, inverse_data[[object@eq_constr]])
      leq_dual <- get_dual_values(solution$ineq_dual, extract_dual_value, inverse_data[[object@neq_constr]])
      eq_dual <- utils::modifyList(eq_dual, leq_dual)
      dual_vars <- eq_dual
    }
  } else {
    primal_vars <- list()
    primal_vars[[as.character(inverse_data[[as.character(object@var_id)]])]] <- NA_real_
    if(!inverse_data[["is_mip"]]) {
      dual_var_ids <- sapply(c(inverse_data[[object@eq_constr]], inverse_data[[object@neq_constr]]), function(constr) { constr@id })
      dual_vars <- as.list(rep(NA_real_, length(dual_var_ids)))
      names(dual_vars) <- dual_var_ids
    }

    if(status == INFEASIBLE)
      opt_val <- Inf
    else if(status == UNBOUNDED)
      opt_val <- -Inf
    else
      opt_val <- NA_real_
  }

  return(Solution(status, opt_val, primal_vars, dual_vars, list()))
})

setMethod("solve_via_data", "GUROBI_CONIC", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
  requireNamespace("gurobi", quietly = TRUE)

  cvar <- data[[C_KEY]]
  b <- data[[B_KEY]]
  A <- data[[A_KEY]]
  dims <- data[[DIMS]]

  n <- length(cvar)

  #Create a new model and add objective term
  model <- list()
  model$obj <- c(cvar, rep(0, sum(unlist(dims@soc))))

  #Add variable types
  vtype <- character(n)
  for(i in seq_along(data[[BOOL_IDX]])){
    vtype[data[[BOOL_IDX]][[i]]] <- 'B' #B for binary
  }
  for(i in seq_along(data[[INT_IDX]])){
    vtype[data[[INT_IDX]][[i]]] <- 'I' #I for integer
  }

  for(i in 1:n) {
    if(vtype[i] == ""){
      vtype[i] <- 'C' #C for continuous
    }
  }
  model$vtype <- vtype #put in variable types
  model$lb <- rep(-Inf, n)
  model$ub <- rep(Inf, n)

  # Add equality constraints: iterate over the rows of A,
  # adding each row into the model.
  model$A <- A
  model$rhs <- b
  model$sense <- c(rep('=', dims@zero), rep('<',  dims@nonpos), rep('=', sum(unlist(dims@soc))))

  total_soc <- sum(unlist(dims@soc))
  current_vars <- n
  current_rows <- dims@zero + dims@nonpos + 1

  # Add SOC variables
  # Sort of strange. A good example of how it works can be seen in
  # https://www.gurobi.com/documentation/8.1/examples/qcp_r.html#subsubsection:qcp.R
  for(i in seq_along(dims@soc)){
    n_soc <- dims@soc[[i]]

    model$vtype <- c(model$vtype, rep('C', n_soc))
    model$lb <- c(model$lb, 0, rep(-Inf, n_soc - 1))
    model$ub <- c(model$ub, rep(Inf, n_soc))
    Asub <- matrix(0, nrow = nrow(A), ncol = n_soc)
    Asub[current_rows:(current_rows + n_soc - 1),] <- diag(rep(1, n_soc))
    model$A <- cbind(model$A, Asub)

    # To create quadratic constraints, first create a 0 square matrix with dimension of
    # the total number of variables (normal + SOC). Then fill the diagonals of the
    # SOC part with the first being negative and the rest being positive
    qc <- list()
    qc$Qc <- matrix(0, nrow = n + total_soc, ncol = n + total_soc)
    qc$Qc[current_vars + 1, current_vars + 1] <- -1
    for(j in 1:(n_soc-1)){
      qc$Qc[current_vars + 1 + j, current_vars + 1 + j] <- 1
    }
    qc$rhs <- 0.0

    model$quadcon[[i]] <- qc

    current_vars <- current_vars + n_soc
    current_rows = current_rows + n_soc

  }

  params <- list()
  params$OutputFlag <- as.numeric(verbose)
  params$QCPDual <- 1 #equivalent to TRUE
  for(i in seq_along(solver_opts)){
    params[[ names(solver_opts)[i] ]] <- solver_opts[i]
  }

  solution <- list()
  tryCatch({
    result <- gurobi::gurobi(model, params)   # Solve.
    solution[["value"]] <- result$objval
    solution[["primal"]] <- result$x

    #Only add duals if it's not a MIP
    if(sum(unlist(data[[BOOL_IDX]])) + sum(unlist(data[[INT_IDX]])) == 0){
      solution[["y"]] <- -append(result$pi, result$qcpi, dims@zero + dims@nonpos)

      if(dims@zero == 0){
        solution[["eq_dual"]] <- c()
        solution[["ineq_dual"]] <- solution[["y"]]
      } else {
        solution[["eq_dual"]] <- solution[["y"]][1:dims@zero]
        solution[["ineq_dual"]] <- solution[["y"]][-(1:dims@zero)]
      }
    }

  }, error = function(e) {   # Error in the solution.
  })

  solution[[SOLVE_TIME]] <- result$runtime
  solution[["status"]] <- status_map(object, result$status)
  solution[["num_iters"]] <- result$baritercount

  # Is there a better way to check if there is a solution?
  # if(solution[["status"]] == SOLVER_ERROR && !is.na(result$x)){
  #  solution[["status"]] <- OPTIMAL_INACCURATE
  # }

  return(solution)

})

MOSEK <- setClass("MOSEK", representation(exp_cone_order = "numeric"),   # Order of exponential cone constraints. Internal only!
                           prototype(exp_cone_order = c(2, 1, 0)), contains = "ConicSolver")

vectorized_lower_tri_to_mat <- function(v, dim) {
  v <- unlist(v)
  rows <- c()
  cols <- c()
  vals <- c()
  running_idx <- 1
  for(j in seq_len(dim)) {
    rows <- c(rows, j + seq_len(dim-j+1) - 1)
    cols <- c(cols, rep(j, dim-j+1))
    vals <- c(vals, v[running_idx:(running_idx + dim - j)])
    running_idx <- running_idx + dim - j + 1
  }
  A <- sparseMatrix(i = rows, j = cols, x = vals, dims = c(dim, dim))
  d <- diag(diag(A))
  A <- A + t(A) - d
  return(A)
}

psd_coeff_offset <- function(problem, c) {
  # Returns an array G and vector h such that the given constraint is
  # equivalent to G*z <=_{PSD} h.
  extractor <- CoeffExtractor(InverseData(problem))
  tmp <- affine(extractor, expr(c))
  A_vec <- tmp[[1]]
  b_vec <- tmp[[2]]
  G <- -A_vec
  h <- b_vec
  dim <- nrow(expr(c))
  return(list(G, h, dim))
}

setMethod("mip_capable", "MOSEK", function(solver) { TRUE })
setMethod("supported_constraints", "MOSEK", function(solver) { c(supported_constraints(ConicSolver()), "SOC", "PSDConstraint") })

setMethod("import_solver", "MOSEK", function(solver) {
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
    return(list(NULL, NULL))
  matrices <- list()
  offsets <- c()
  lengths <- c()
  ids <- c()

  for(con in constraints) {
    coeff_offs <- reduction_format_constr(object, problem, con, exp_cone_order)
    coeff <- coeff_offs[[1]]
    offset <- coeff_offs[[2]]
    matrices <- c(matrices, list(coeff))
    offsets <- c(offsets, offset)
    lengths <- c(lengths, prod(dim(offset)))
    ids <- c(ids, id(con))
  }
  coeff <- Matrix(do.call(rbind, matrices), sparse = TRUE)
  return(list(coeff, offsets, lengths, ids))
})

setMethod("perform", signature(object = "MOSEK", problem = "Problem"), function(object, problem) {
  data <- list()
  inv_data <- list(suc_slacks = list(), y_slacks = list(), snx_slacks = list(), psd_dims = list())
  inv_data[[object@var_id]] <- id(variables(problem)[[1]])

  # Get integrality constraint information.
  var <- variables(problem)[[1]]
  data[[BOOL_IDX]] <- sapply(var@boolean_idx, function(t) { as.integer(t[1]) })
  data[[INT_IDX]] <- sapply(var@integer_idx, function(t) { as.integer(t[1]) })
  inv_data$integer_variables <- length(data[[BOOL_IDX]]) + length(data[[INT_IDX]]) > 0

  # Parse the coefficient vector from the objective.
  coeff_offs <- ConicSolver.get_coeff_offset(problem@objective@args[[1]])
  c <- coeff_offs[[1]]
  constant <- coeff_offs[[2]]
  data[[C_KEY]] <- as.vector(c)
  inv_data$n0 <- length(data[[C_KEY]])
  data[[OBJ_OFFSET]] <- constant[1]
  data[[DIMS]] <- list()
  data[[DIMS]][[SOC_DIM]] <- list()
  data[[DIMS]][[EXP_DIM]] <- list()
  data[[DIMS]][[PSD_DIM]] <- list()
  data[[DIMS]][[LEQ_DIM]] <- 0
  data[[DIMS]][[EQ_DIM]] <- 0
  inv_data[[OBJ_OFFSET]] <- constant[1]
  Gs <- list()
  hs <- list()

  if(length(problem@constraints) == 0) {
    ##data[[G_KEY]] <- Matrix(nrow = 0, ncol = 0, sparse = TRUE)
    ## Ensure G's dimensions match that of c.
    data[[G_KEY]] <- Matrix(nrow = 0, ncol = length(c), sparse = TRUE)
    data[[H_KEY]] <- matrix(nrow = 0, ncol = 0)
    inv_data$is_LP <- TRUE
    return(list(object, data, inv_data))
  }

  # Linear inequalities.
  leq_constr <- problem@constraints[sapply(problem@constraints, function(ci) { class(ci) == "NonPosConstraint" })]
  if(length(leq_constr) > 0) {
    blform <- block_format(object, problem, leq_constr)   # G, h : G*z <= h.
    G <- blform[[1]]
    h <- blform[[2]]
    lengths <- blform[[3]]
    ids <- blform[[4]]
    inv_data$suc_slacks <- c(inv_data$suc_slacks, lapply(1:length(lengths), function(k) { c(ids[k], lengths[k]) }))
    data[[DIMS]][[LEQ_DIM]] <- sum(lengths)
    Gs <- c(Gs, G)
    hs <- c(hs, h)
  }

  # Linear equations.
  eq_constr <- problem@constraints[sapply(problem@constraints, function(ci) { class(ci) == "ZeroConstraint" })]
  if(length(eq_constr) > 0) {
    blform <- block_format(object, problem, eq_constr)   # G, h : G*z == h.
    G <- blform[[1]]
    h <- blform[[2]]
    lengths <- blform[[3]]
    ids <- blform[[4]]
    inv_data$y_slacks <- c(inv_data$y_slacks, lapply(1:length(lengths), function(k) { c(ids[k], lengths[k]) }))
    data[[DIMS]][[EQ_DIM]] <- sum(lengths)
    Gs <- c(Gs, G)
    hs <- c(hs, h)
  }

  # Second order cone.
  soc_constr <- problem@constraints[sapply(problem@constraints, function(ci) { class(ci) == "SOC" })]
  data[[DIMS]][[SOC_DIM]] <- list()
  for(ci in soc_constr)
    data[[DIMS]][[SOC_DIM]] <- c(data[[DIMS]][[SOC_DIM]], cone_sizes(ci))
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
    blform <- block_format(object, problem, exp_constr, object@exp_cone_order)
    G <- blform[[1]]
    h <- blform[[2]]
    lengths <- blform[[3]]
    ids <- blform[[4]]
    data[[DIMS]][[EXP_DIM]] <- lengths
    Gs <- c(Gs, G)
    hs <- c(hs, h)
  }

  # PSD constraints.
  psd_constr <- problem@constraints[sapply(problem@constraints, function(ci) { class(ci) == "PSDConstraint" })]
  if(length(psd_constr) > 0) {
    data[[DIMS]][[PSD_DIM]] <- list()
    for(c in psd_constr) {
      coeff_offs <- psd_coeff_offset(problem, c)
      G_vec <- coeff_offs[[1]]
      h_vec <- coeff_offs[[2]]
      dim <- coeff_offs[[3]]
      inv_data$psd_dims <- c(inv_data$psd_dims, list(list(id(c), dim)))
      data[[DIMS]][[PSD_DIM]] <- c(data[[DIMS]][[PSD_DIM]], list(dim))
      Gs <- c(Gs, G_vec)
      hs <- c(hs, h_vec)
    }
  }

  if(length(Gs) == 0)
    ## data[[G_KEY]] <- Matrix(nrow = 0, ncol = 0, sparse = TRUE)
    ## G is already sparse
    data[[G_KEY]] <- G
  else
    data[[G_KEY]] <- Matrix(do.call(rbind, Gs), sparse = TRUE)
  if(length(hs) == 0)
    data[[H_KEY]] <- matrix(nrow = 0, ncol = 0)
  else
    data[[H_KEY]] <- Matrix(do.call(cbind, hs), sparse = TRUE)
  inv_data$is_LP <- (length(psd_constr) + length(exp_constr) + length(soc_constr)) == 0
  return(list(object, data, inv_data))
})

# TODO: Finish MOSEK class implementation.
setMethod("solve_via_data", "MOSEK", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
    ##requireNamespace("Rmosek", quietly = TRUE)
  #env <- Rmosek::Env() remove these as Rmosek doesn't need environments
  #task <- env.Task(0,0)
  #instead defines prob
  prob <- list(sense="min")

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
  if(length(data[[C_KEY]]) == 0) {
    res <- list()
    res[[STATUS]] <- OPTIMAL
    res[[PRIMAL]] <- list()
    res[[VALUE]] <- data[[OFFSET]]
    res[[EQ_DUAL]] <- list()
    res[[INEQ_DUAL]] <- list()
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

  c <- data[[C_KEY]]
  G <- data[[G_KEY]]
  h <- data[[H_KEY]]
  dims <- data[[DIMS]]
  n0 <- length(c)
  n <- n0 + sum(unlist(dims[[SOC_DIM]]), na.rm = TRUE) + sum(unlist(dims[[EXP_DIM]]), na.rm = TRUE) # unlisted dims to make sure sum function works and na.rm to handle empty lists
  psd_total_dims <- sum(unlist(dims[[PSD_DIM]])^2, na.rm = TRUE)
  m <- length(h)
  num_bool <- length(data[[BOOL_IDX]])
  num_int <- length(data[[INT_IDX]])

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


  #task.appendvars(n), no task for Rmosek, but declares the number of variables in the model. Need to expand prob$c as well to match this dimension
  #task.putvarboundlist(1:n, rep(mosek.boundkey.fr, n), matrix(0, nrow = n, ncol = 1), matrix(0, nrow = n, ncol = 1))
  #kind of confused why x's are all 0's, but that's what's in python code
  #prob$bx <- rbind( blx = rep(0,n),
  #                 bux = rep(0,n))
  prob$bx <- rbind(blx = rep(-Inf, n),
                   bux = rep(Inf, n))

  #Initialize the cone. Not 100% sure about this bit
  NUMCONES <- length(dims[[SOC_DIM]]) + floor(sum(unlist(dims[[EXP_DIM]]), na.rm = TRUE)/3)
  prob$cones <- matrix(list(), nrow = 2, ncol = NUMCONES)

  if(psd_total_dims > 0)
    prob$bardim <- unlist(dims[[PSD_DIM]])
  running_idx <- n0
  for(i in seq_along(unlist(dims[[SOC_DIM]]))) {
    prob$cones[,i] <- list("QUAD", as.numeric((running_idx + 1):(running_idx + unlist(dims[[SOC_DIM]])[[i]]))) # latter term is size_cone
    running_idx <- running_idx + unlist(dims[[SOC_DIM]])[[i]]
  }
  if(floor(sum(unlist(dims[[EXP_DIM]]), na.rm = TRUE)/3) != 0){ # check this, feels sketchy
    for(k in 1:(floor(sum(unlist(dims[[EXP_DIM]]), na.rm = TRUE)/3)+1) ) {
      prob$cones[,(length(dims[[SOC_DIM]])+k)] <- list("PEXP", as.numeric((running_idx+1):(running_idx + 3)) )
      running_idx <- running_idx + 3
    }
  }
  if(num_bool + num_int > 0) {
    if(num_bool > 0) {
      prob$intsub <- unlist(data[[BOOL_IDX]])
      #since the variable constraints are already declared, we are resetting them so they can only be 0 or 1
      prob$bx[, unlist(data[[BOOL_IDX]])] <- rbind( rep(0, length(unlist(data[[BOOL_IDX]]))), rep(1, length(unlist(data[[BOOL_IDX]]))) )

    }
    if(num_int > 0)
      prob$intsub <- unlist(data[[INT_IDX]])
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

  # task.appendcons(m) is equivalent to prob$bc


  ##G_sparse <- as(as.matrix(G), "sparseMatrix")
  ##G is already sparse
  G_sparse  <- G
  G_sum <- summary(G_sparse)
  row <- G_sum$i
  col <- G_sum$j
  vals <- G_sum$x

  total_soc_exp_slacks <- sum(unlist(dims[[SOC_DIM]]), na.rm = TRUE) + sum(unlist(dims[[EXP_DIM]]), na.rm = TRUE)

  # initializing A matrix
  if(nrow(G_sparse) == 0 || (ncol(G_sparse) + total_soc_exp_slacks) == 0)
    ## prob$A <- sparseMatrix(i = c(), j = c(), dims = c(0, 0))
    ## G is already sparse
    prob$A  <- G
  else {
    prob$A <- sparseMatrix(i = rep(1:nrow(G_sparse), ncol(G_sparse) + total_soc_exp_slacks),
                           j = rep(1:(ncol(G_sparse) + total_soc_exp_slacks), nrow(G_sparse)),
                           x = rep(0, nrow(G_sparse)*(ncol(G_sparse) + total_soc_exp_slacks)))

    # this is a bit hacky, probably should fix later. Filling out part of the A matrix from G
    # Equivalent to task.putaijlist(as.list(row), as.list(col), as.list(vals))
    A_holder <- sparseMatrix(row, col, x = vals)
    prob$A[1:nrow(A_holder), 1:ncol(A_holder)] <- A_holder
  }

  if(total_soc_exp_slacks > 0) {
    i <- unlist(dims[[LEQ_DIM]]) + unlist(dims[[EQ_DIM]])   # Constraint index in (1, ..., m)
    j <- length(c)   # Index of the first slack variable in the block vector "x".
    rows <- (i:(i + total_soc_exp_slacks-1))+1
    cols <- (j:(j + total_soc_exp_slacks-1))+1
    #task.putaijlist(rows, cols, rep(1, total_soc_exp_slacks))
    for(iter in 1:length(rows))
      prob$A[rows[iter],cols[iter]] <- 1
  }

  # Constraint index: start of LMIs.
  i <- dims[[LEQ_DIM]] + dims[[EQ_DIM]] + total_soc_exp_slacks + 1
  dim_exist_PSD <- length(dims[[PSD_DIM]]) #indicates whether or not we have any LMIs

  if(dim_exist_PSD > 0){
    #A bit hacky here too, specifying the lower triangular part of symmetric coefficient matrix barA
    barAi <- c() #Specifies row index of block matrix
    barAj <- c() #Specifies column index of block matrix
    barAk <- c() #Specifies row index within the block matrix specified above
    barAl <- c() #Specifies column index within the block matrix specified above
    barAv <- c() #Values for all the matrices

    for(j in 1:length(dims[[PSD_DIM]])) {   #For each PSD matrix
      for(row_idx in 1:dims[[PSD_DIM]][[j]]) {
        for(col_idx in 1:dims[[PSD_DIM]][[j]]) {
          val <- ifelse(row_idx == col_idx, 1, 0.5)
          row <- max(row_idx, col_idx)
          col <- min(row_idx, col_idx)
          #mat <- task.appendsparsesymmat(dim, list(row), list(col), list(val))
          #task.putbaraij(i, j, list(mat), list(1.0))

          barAi <- c(barAi, i)
          barAj <- c(barAj, j) #NEED TO CHECK. Multiple PSD_DIM example?
          barAk <- c(barAk, row)
          barAl <- c(barAl, col)
          barAv <- c(barAv, val)

          i <- i + 1 #for each symmetric matrix
        }
      }
    }

    #Attaching. Does mosek automatically check the symmetric matrix dimensions?

    prob$barA$i <- barAi
    prob$barA$j <- barAj
    prob$barA$k <- barAk
    prob$barA$l <- barAl
    prob$barA$v <- barAv
  }

  num_eq <- length(h) - dims[[LEQ_DIM]]

  #CVXPY has the first dims[[LEQ_DIM]] variables as upper bounded
  #type_constraint <- rep(mosek.boundkey.up, dims[[LEQ_DIM]]) + rep(mosek.boundkey.fx, num_eq)
  #task.putconboundlist(1:m, type_constraint, h, h), equivalent to prob$bc
  hl_holder <- as.numeric(h)
  hu_holder <- as.numeric(h)

  #upper constraints for the LEQ_DIM number of variables, so set lower bound to -Inf
  hl_holder[seq_len(dims[[LEQ_DIM]])] <- rep(-Inf, dims[[LEQ_DIM]])

  prob$bc <- rbind(blc = hl_holder,
                   buc = hu_holder)

  # Define the objective and optimize the MOSEK task.
  #task.putclist(1:length(c), c)
  #initialize coefficients of objective with the same number of variables declared (dim of x)
  c_holder <- rep(0, n)
  c_holder[1:length(c)] <- c

  prob$c <- c_holder
  #task.putobjsense(mosek.objsense.minimize) Rmosek does this at the beginning of the problem instead of the end like python

  #if(!is.na(save_file)) Don't think there's a save equivalent in R
  #  task.writedata(save_file)
  #task.optimize

  if(verbose){
    r <- Rmosek::mosek(prob, list(verbose = 1, soldetail = 3, getinfo = TRUE))
    ##sol <- r$sol
  } else {
    r <- Rmosek::mosek(prob, list(soldetail = 3, getinfo = TRUE))
    ##sol <- r$sol
  }
  r
##  return(sol)
})

setMethod("invert", "MOSEK", function(object, solution, inverse_data) {
  ## REMOVE LATER
  ## results  <- solution
  ##    has_attr <- !is.null(mosek.solsta$near_optimal)
  ## We ignore MOSEK 8.1 and below.
  status_map <- function(status) {
      status  <- tolower(status)
      if(status %in% c("optimal", "integer_optimal"))
          return(OPTIMAL)
      ##        else if(status %in% c("prim_feas", "near_optimal", "near_integer_optimal"))
      ##            return(OPTIMAL_INACCURATE)
      else if(status == "prim_infeas_cer" || status == "primal_infeasible_cer") { #Documentation says it's this, but docs also say it spits out dual_infeas_cer, which is wrong
        #check later
          if(!is.null(attributes(status))) #check if status has any attributes, hasattr in python
              return(INFEASIBLE)
          else
              return(INFEASIBLE)
      } else if(status == "dual_infeasible_cer") {
          if(!is.null(attributes(status)))
              return(UNBOUNDED_INACCURATE)
          else
              return(UNBOUNDED)
      }  else
          return(SOLVER_ERROR)
  }

  ##env <- results$env
  ##task <- results$task
  ## Naras: FIX solver_opts
  solver_opts <- solution$solver_options

  if(inverse_data$integer_variables)
      sol <- solution$sol$int
  else if(!is.null(solver_opts$bfs) && solver_opts$bfs && inverse_data$is_LP)
      sol <- solution$sol$bas   # The basic feasible solution.
  else
      sol <- solution$sol$itr   # The solution found via interior point method.

  problem_status <- sol$prosta
  solution_status <- sol$solsta

  if(is.na(solution$response$code))
    status <- SOLVER_ERROR
  else
    status <- status_map(solution_status)

  ## For integer problems, problem status determines infeasibility (no solution).
  ##  if(sol == mosek.soltype.itg && problem_status == mosek.prosta.prim_infeas)
  ## Using reference https://docs.mosek.com/9.0/rmosek/accessing-solution.html
  if(inverse_data$integer_variables && (problem_status == "MSK_PRO_STA_PRIM_INFEAS" || problem_status == "PRIMAL_INFEASIBLE"))
      status <- INFEASIBLE

  if(status %in% SOLUTION_PRESENT) {
                                      # Get objective value.
      opt_val <- sol$pobjval + inverse_data[[OBJ_OFFSET]]
                                      # Recover the CVXR standard form primal variable.
      ## z <- rep(0, inverse_data$n0)
      ## task.getxxslice(sol, 0, length(z), z)
      primal_vars <- list()
      primal_vars[[as.character(inverse_data[[object@var_id]])]] <- sol$xx
      ## Recover the CVXR standard form dual variables.
      ## if(sol == mosek.soltype.itn)
      if (inverse_data$integer_variables) {
        dual_var_ids <- sapply(c(inverse_data$suc_slacks, inverse_data$y_slacks, inverse_data$snx_slacks, inverse_data$psd_dims), function(slack) { slack[[1L]] })
        dual_vars <- as.list(rep(NA_real_, length(dual_var_ids)))
        names(dual_vars) <- dual_var_ids
      } else
        dual_vars <- MOSEK.recover_dual_variables(task, sol, inverse_data)

  } else {
      if(status == INFEASIBLE)
          opt_val <- Inf
      else if(status == UNBOUNDED)
          opt_val <- -Inf
      else
          opt_val <- NA_real_
      vid <-
      primal_vars <- list()
      primal_vars[[as.character(inverse_data[[object@var_id]])]] <- NA_real_
      dual_var_ids <- sapply(c(inverse_data$suc_slacks, inverse_data$y_slacks, inverse_data$snx_slacks, inverse_data$psd_dims), function(slack) { slack[[1L]] })
      dual_vars <- as.list(rep(NA_real_, length(dual_var_ids)))
      names(dual_vars) <- dual_var_ids
  }

  ## Store computation time.
  attr <- list()
  attr[[SOLVE_TIME]] <- solution$dinfo$OPTIMIZER_TIME

  ## Delete the MOSEK Task and Environment
  ##task.__exit__(NA, NA, NA)
  ##env.__exit__(NA, NA, NA)

  return(Solution(status, opt_val, primal_vars, dual_vars, attr))
})

MOSEK.recover_dual_variables <- function(task, sol, inverse_data) {
  dual_vars <- list()

  ## Dual variables for the inequality constraints.
  suc_len <- ifelse(length(inverse_data$suc_slacks) == 0, 0, sum(sapply(inverse_data$suc_slacks, function(val) { val[[2]] })))
  if(suc_len > 0) {
      ## suc <- rep(0, suc_len)
      ## task.getsucslice(sol, 0, suc_len, suc)
      dual_vars <- utils::modifyList(dual_vars, MOSEK.parse_dual_vars(sol$suc[seq_len(suc_len)], inverse_data$suc_slacks))
  }

  ## Dual variables for the original equality constraints.
  y_len <- ifelse(length(inverse_data$y_slacks) == 0, 0, sum(sapply(inverse_data$y_slacks, function(val) { val[[2]] })))
  if(y_len > 0) {
      ##y <- rep(0, y_len)
      ## task.getyslice(sol, suc_len, suc_len + y_len, y)
      dual_vars <- utils::modifyList(dual_vars, MOSEK.parse_dual_vars(sol$suc[seq.int(suc_len, length.out = y_len)], inverse_data$y_slacks))
  }

  ## Dual variables for SOC and EXP constraints.
  snx_len <- ifelse(length(inverse_data$snx_slacks) == 0, 0, sum(sapply(inverse_data$snx_slacks, function(val) { val[[2]] })))
  if(snx_len > 0) {
      ##snx <- matrix(0, nrow = snx_len, ncol = 1)
      ##task.getsnxslice(sol, inverse_data$n0, inverse_data$n0 + snx_len, snx)
      dual_vars <- utils::modifyList(dual_vars, MOSEK.parse_dual_vars(sol$snx, inverse_data$snx_slacks))
  }

  ## Dual variables for PSD constraints.
  for(psd_info in inverse_data$psd_dims) {
    id <- as.character(psd_info[[1L]])
    dim <- psd_info[[2L]]
    ##sj <- rep(0, dim*floor((dim + 1)/2))
    ##task.getbars(sol, j, sj)
    dual_vars[[id]] <- vectorized_lower_tri_to_mat(sol$bars[[1L]], dim)
  }

  return(dual_vars)
}

MOSEK.parse_dual_vars <- function(dual_var, constr_id_to_constr_dim) {
  dual_vars <- list()
  running_idx <- 1
  for(val in constr_id_to_constr_dim) {
    id <- as.character(val[[1]])
    dim <- val[[2]]
    ## if(dim == 1)
    ##   dual_vars[id] <- dual_vars[running_idx]   # a scalar.
    ## else
    ##   dual_vars[id] <- as.matrix(dual_vars[running_idx:(running_idx + dim)])
    dual_vars[[id]] <- dual_var[seq.int(running_idx, length.out = dim)]
    running_idx <- running_idx + dim
  }
  return(dual_vars)
}

MOSEK._handle_mosek_params <- function(task, params) {
  if(is.na(params))
    return()

  ##requireNamespace("Rmosek", quietly = TRUE)

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
SCS.dims_to_solver_dict <- function(cone_dims) {
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
  entries <- floor(rows * (cols + 1)/2)

  row_arr <- seq_len(entries)

  col_arr <- matrix(1:(rows*cols), nrow = rows, ncol = cols)
  col_arr <- col_arr[lower.tri(col_arr, diag = TRUE)]

  val_arr <- matrix(0, nrow = rows, ncol = cols)
  val_arr[lower.tri(val_arr, diag = TRUE)] <- sqrt(2)
  diag(val_arr) <- 1
  val_arr <- as.vector(val_arr)
  val_arr <- val_arr[val_arr != 0]

  coeff <- Constant(sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = c(entries, rows*cols)))
  vectorized_matrix <- reshape_expr(matrix, c(rows*cols, 1))
  return(coeff %*% vectorized_matrix)
}

tri_to_full <- function(lower_tri, n) {
  # Expands n*floor((n+1)/2) lower triangular to full matrix.
  # Scales off-diagonal by 1/sqrt(2), as per the SCS specification.
  full <- matrix(0, nrow = n, ncol = n)
  full[upper.tri(full, diag = TRUE)] <- lower_tri
  full[lower.tri(full, diag = TRUE)] <- lower_tri

  unscaled_diag <- diag(full)
  full <- full/sqrt(2)
  diag(full) <- unscaled_diag

  matrix(full, nrow = n*n, byrow = FALSE)
}

SuperSCS <- setClass("SuperSCS", contains = "SCS")
SuperSCS.default_settings <- function(object) {
  list(use_indirect = FALSE, eps = 1e-8, max_iters = 10000)
}

setMethod("name", "SuperSCS", function(x) { SUPER_SCS_NAME })
setMethod("import_solver", "SuperSCS", function(solver) {
  stop("Unimplemented: SuperSCS is currently unavailable in R.")
})

setMethod("solve_via_data", "SuperSCS", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
  args <- list(A = data[[A_KEY]], b = data[[B_KEY]], c = data[[C_KEY]])
  if(warm_start && !is.null(solver_cache) && length(solver_cache) > 0 && name(object) %in% names(solver_cache)) {
    args$x <- solver_cache[[name(object)]]$x
    args$y <- solver_cache[[name(object)]]$y
    args$s <- solver_cache[[name(object)]]$s
  }
  cones <- SCS.dims_to_solver_dict(data[[ConicSolver()@dims]])

  # Settings.
  user_opts <- names(solver_opts)
  for(k in names(SuperSCS.default_settings)) {
    if(!k %in% user_opts)
      solver_opts[[k]] <- SuperSCS.default_settings[[k]]
  }
  results <- SuperSCS::solve(args, cones, verbose = verbose, solver_opts)
  if(!is.null(solver_cache) && length(solver_cache) > 0)
    solver_cache[[name(object)]] <- results
  return(results)
})

# XPRESS <- setClass("XPRESS", contains = "SCS")
#
# # Solver capabilities.
# setMethod("mip_capable", "XPRESS", function(solver) { TRUE })
# setMethod("supported_constraints", "XPRESS", function(solver) { c(supported_constraints(ConicSolver()), "SOC") })
#
# # Map of XPRESS status to CVXR status.
# setMethod("status_map", "XPRESS", function(solver, status) {
#   if(status == 2)
#     return(OPTIMAL)
#   else if(status == 3)
#     return(INFEASIBLE)
#   else if(status == 5)
#     return(UNBOUNDED)
#   else if(status %in% c(4, 6, 7, 8, 10, 11, 12, 13))
#     return(SOLVER_ERROR)
#   else if(status == 9)   # TODO: Could be anything. Means time expired.
#     return(OPTIMAL_INACCURATE)
#   else
#     stop("XPRESS status unrecognized: ", status)
# })
#
# setMethod("name", "XPRESS", function(x) { XPRESS_NAME })
# setMethod("import_solver", "XPRESS", function(solver) {
#   stop("Unimplemented: XPRESS solver unavailable in R.")
# })
#
# setMethod("accepts", signature(object = "XPRESS", problem = "Problem"), function(object, problem) {
#   # TODO: Check if the matrix is stuffed.
#   if(!is_affine(problem@objective@args[[1]]))
#     return(FALSE)
#   for(constr in problem@constraints) {
#     if(!class(constr) %in% supported_constraints(object))
#       return(FALSE)
#     for(arg in constr@args) {
#       if(!is_affine(arg))
#         return(FALSE)
#     }
#   }
#   return(TRUE)
# })
#
# setMethod("perform", signature(object = "XPRESS", problem = "Problem"), function(object, problem) {
#   tmp <- callNextMethod(object, problem)
#   data <- tmp[[1]]
#   inv_data <- tmp[[2]]
#   variables <- variables(problem)[[1]]
#   data[[BOOL_IDX]] <- lapply(variables@boolean_idx, function(t) { t[1] })
#   data[[INT_IDX]] <- lapply(variables@integer_idx, function(t) { t[1] })
#   inv_data$is_mip <- length(data[[BOOL_IDX]]) > 0 || length(data[[INT_IDX]]) > 0
#   return(list(object, data, inv_data))
# })
#
# setMethod("invert", signature(object = "XPRESS", solution = "list", inverse_data = "list"), function(object, solution, inverse_data) {
#   status <- solution[[STATUS]]
#
#   if(status %in% SOLUTION_PRESENT) {
#     opt_val <- solution[[VALUE]]
#     primal_vars <- list()
#     primal_vars[[inverse_data[[object@var_id]]]] <- solution$primal
#     if(!inverse_data@is_mip)
#       dual_vars <- get_dual_values(solution[[EQ_DUAL]], extract_dual_value, inverse_data[[EQ_CONSTR]])
#   } else {
#     primal_vars <- list()
#     primal_vars[[inverse_data[[object@var_id]]]] <- NA_real_
#     if(!inverse_data@is_mip) {
#       dual_var_ids <- sapply(inverse_data[[EQ_CONSTR]], function(constr) { constr@id })
#       dual_vars <- as.list(rep(NA_real_, length(dual_var_ids)))
#       names(dual_vars) <- dual_var_ids
#     }
#
#     if(status == INFEASIBLE)
#       opt_val <- Inf
#     else if(status == UNBOUNDED)
#       opt_val <- -Inf
#     else
#       opt_val <- NA
#   }
#
#   other <- list()
#   other[[XPRESS_IIS]] <- solution[[XPRESS_IIS]]
#   other[[XPRESS_TROW]] <- solution[[XPRESS_TROW]]
#   return(Solution(status, opt_val, primal_vars, dual_vars, other))
# })
#
# setMethod("solve_via_data", "XPRESS", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
#   solver <- XPRESS_OLD()
#   solver_opts[[BOOL_IDX]] <- data[[BOOL_IDX]]
#   solver_opts[[INT_IDX]] <- data[[INT_IDX]]
#   prob_data <- list()
#   prob_data[[name(object)]] <- ProblemData()
#   solve(solver, data$objective, data$constraints, prob_data, warm_start, verbose, solver_opts)
# })
