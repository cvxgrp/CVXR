# QP interface for the CPLEX solver.
CPLEX <- setClass("CPLEX", contains = "QpSolver")

setMethod("mip_capable", "CPLEX", function(object) { TRUE })

# Map of CPLEX status to CVXR status.
# TODO: Add more!
setMethod("status_map", "OSQP", function(solver, status, default = NULL) {
  if(status %in% c(1, 101))
    OPTIMAL
  else if(status %in% c(3, 22, 4, 103))
    INFEASIBLE
  else if(status %in% c(2, 21, 118))
    UNBOUNDED
  else if(status %in% c(10, 107))
    USER_LIMIT
  else if(!is.null(default)) {
    warning("CPLEX status unrecognized: ", status)
    return(default)
  } else
    stop("CPLEX status unrecognized: ", status)
})

setMethod("name", "CPLEX", function(x) { CPLEX_NAME })
setMethod("import_solver", "CPLEX", function(object) { requireNamespace(Rcplex, quietly = TRUE) })
setMethod("invert", "CPLEX", function(object, results, inverse_data) {
  model <- results$model
  attr <- list()
  if("cputime" %in% names(results))
    attr[SOLVE_TIME] <- results$cputime
  attr[NUM_ITERS] <- as.integer(get_num_barrier_iterations(model@solution@progress))
  
  status <- status_map(object, get_status(model@solution), SOLVER_ERROR)
  
  if(status %in% SOLUTION_PRESENT) {
    # Get objective value.
    opt_val <- get_objective_value(model@solution)
    
    # Get solution.
    x <- as.matrix(get_values(model@solution))
    primal_vars <- list()
    primal_vars[names(inverse_data@id_map)[1]] <- x
    
    # Only add duals if not a MIP.
    dual_vars <- NA
    if(!inverse_data@is_mip) {
      y <- -as.matrix(get_dual_values(model@solution))
      dual_vars <- get_dual_values(y, extract_dual_value, inverse_data@sorted_constraints)
    }
  } else {
    primal_vars <- NA
    dual_vars <- NA
    opt_val <- Inf
    if(status == UNBOUNDED)
      opt_val <- -Inf
  }
  
  return(Solution(status, opt_val, primal_vars, dual_vars, attr))
})

setMethod("solve_via_data", "CPLEX", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  requireNamespace(Rcplex, quietly = TRUE)
  # TODO: Finish this.
})

# QP interface for the OSQP solver.
OSQP <- setClass("OSQP", contains = "QpSolver")

# Map of OSQP status to CVXPY status.
setMethod("status_map", "OSQP", function(solver, status, default = NULL) {
  if(status == 1)
    OPTIMAL
  else if(status == 2)
    OPTIMAL_INACCURATE
  else if(status == -2)   # Maxiter reached.
    SOLVER_ERROR
  else if(status == -3)
    INFEASIBLE
  else if(status == 3)
    INFEASIBLE_INACCURATE
  else if(status == -4)
    UNBOUNDED
  else if(status == 4)
    UNBOUNDED_INACCURATE
  else if(status == -5)   # Interrupted by user.
    SOLVER_ERROR
  else if(status == -10)  # Unsolved.
    SOLVER_ERROR
  else if(!is.null(default)) {
    warning("OSQP status unrecognized: ", status)
    return(default)
  } else
    stop("OSQP status unrecognized: ", status)
})

setMethod("name", "OSQP", function(x) { OSQP_NAME })
setMethod("import_solver", "OSQP", function(object) {
  requireNamespace(osqp, quietly = TRUE)
})

setMethod("invert", signature(object = "OSQP", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  attr <- list()
  attr[SOLVE_TIME] <- solution@info@run_time
  
  # Map OSQP statuses back to CVXR statuses.
  status <- status_map(object, solution@info@status_val, SOLVER_ERROR)
  
  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution@info@obj_val
    primal_vars <- list()
    primal_vars[names(inverse_data@id_map)] <- solution@x
    dual_vars <- get_dual_values(solution@y, extract_dual_value, inverse_data@sorted_constraints)
    attr[NUM_ITERS] <- solution@info@iter
  } else {
    primal_vars <- NA
    dual_vars <- NA
    opt_val <- Inf
    if(status == UNBOUNDED)
      opt_val <- -Inf
  }
  
  return(Solution(status, opt_val, primal_vars, dual_vars, attr))
})

setMethod("solve_via_data", "OSQP", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  requireNamespace(osqp, quietly = TRUE)
  P <- data[P_KEY]
  q <- data[Q_KEY]
  A <- Matrix(do.call(rbind, list(data[A_KEY], data[F_KEY])), sparse = TRUE)
  data$full_A <- A
  uA <- c(data[B_KEY], data[G_KEY])
  data$u <- uA
  lA <- c(data[B_KEY], rep(-Inf, length(data[G_KEY])))
  data$l <- lA
  
  if(!is.na(solver_cache) && name(object) %in% solver_cache) {
    # Use cached data.
    cache <- solver_cache[name(object)]
    solver <- cache[[1]]
    old_data <- cache[[2]]
    results <- cache[[3]]
    
    same_pattern <- all(dim(P) == dim(old_data[P_KEY])) && all(P@p == old_data[P_KEY]@p) &&
                    all(dim(A) == dim(old_data$full_A)) && all(A@p == old_data$full_A@p)
  } else
    same_pattern <- FALSE
  
  # TODO: Check syntax for setting up R's OSQP solver.
  # If sparsity pattern differs, need to do setup.
  if(warm_start && same_pattern) {
    new_args <- list()
    for(key in c("q", "l", "u")) {
      if(any(data[key] != old_data[key]))
        new_args[key] <- data[key]
    }
    factorizing <- FALSE
    
    if(any(P@i != old_data[P_KEY]@i) || any(P@j != old_data[P_KEY]@j)) {
      new_args$Px_idx <- cbind(P@i, P@j)
      factorizing <- TRUE
    }
    if(any(P@x != old_data[P_KEY]@x)) {
      new_args$Px <- P@x
      factorizing <- TRUE
    }
    if(any(A@i != old_data$full_A@i) || any(A@j != old_data$full_A@j)) {
      new_args$Ax_idx <- cbind(A@i, A@j)
      factorizing <- TRUE
    }
    if(any(A@x != old_data$full_A@x)) {
      new_args$Ax <- A@x
      factorizing <- TRUE
    }
    if(length(new_args) > 0)
      update(solver, new_args)
    
    # Map OSQP statuses back to CVXR statuses.
    status <- status_map(object, results@info@status_val, SOLVER_ERROR)
    if(status == OPTIMAL)
      warm_start(solver, results@x, results@y)
    
    # Polish if factorizing.
    if(is.null(solver_opts$polish))
      solver_opts$polish <- factorizing
    update_settings(solver, verbose = verbose, solver_opts)
  } else {
    # Initialize and solve problem.
    if(is.null(solver_opts$polish))
      solver_opts$polish <- TRUE
    solver <- osqp::OSQP()
    setup(solver, P, q, A, lA, uA, verbose = verbose, solver_opts)
  }
  
  results <- solve(solver)
  
  if(!is.na(solver_cache))
    solver_cache[name(object)] <- list(solver, data, results)
  return(results)
})

# QPSolver requires objectives to be stuffed in the following way.
is_stuffed_qp_objective <- function(objective) {
  expr <- objective@expr
  return(class(expr) == "AddExpression" && length(expr@args) == 2 && class(expr@args[[1]]) == "QuadForm" && class(expr@args[[2]]) == "MulExpression" && is_affine(expr@args[[2]]))
}

# A QP solver interface.
QpSolver <- setClass("QpSolver", contains = "ReductionSolver")

setMethod("accepts", signature(object = "QpSolver", problem = "Problem"), function(object, problem) {
  return(class(problem@objective) == "Minimize" && is_stuffed_qp_objective(problem@objective) && are_args_affine(problem@constraints) &&
         all(sapply(problem@constraints, function(c) { class(c) == "Zero" || class(c) == "NonPos" })))
})

setMethod("apply", signature(object = "QpSolver", problem = "Problem"), function(object, problem) {
  # Construct QP problem data stored in a dictionary.
  # The QP has the following form
  #    minimize 1/2 x' P x + q' x
  #    subject to A x = b
  #               F x <= g
  inverse_data <- InverseData(problem)
  
  obj <- problem@objective
  # quadratic part of objective is x.T * P * x, but solvers expect 0.5*x.T * P * x.
  P <- 2*value(obj@expr@args[[1]]@args[[2]])
  q <- as.vector(value(obj@expr@args[[2]]@args[[1]]))
  
  # Get number of variables.
  n <- problem@size_metrics@num_scalar_variables
  
  # TODO: This dependence on ConicSolver is hacky; something should change here.
  eq_cons <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "Zero" })]
  if(length(eq_cons) > 0) {
    eq_coeffs <- list(list(), c())
    for(con in eq_cons) {
      coeff_offset <- ConicSolver.get_coeff_offset(con@expr)
      eq_coeffs[[1]] <- c(eq_coeffs[[1]], coeff_offset[[1]])
      eq_coeffs[[2]] <- c(eq_coeffs[[2]], coeff_offset[[2]])
    }
    A <- Matrix(do.call(rbind, eq_coeffs[[1]]), byrow = TRUE, sparse = TRUE)
    b <- -eq_coeffs[[2]]
  } else {
    A <- Matrix(nrow = 0, ncol = n, byrow = TRUE, sparse = TRUE)
    b <- -matrix(nrow = 0, ncol = 0)
  }
  
  ineq_cons <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "NonPos" })]
  if(length(ineq_cons) > 0) {
    ineq_coeffs <- list(list(), c())
    for(con in ineq_cons) {
      coeff_offset <- ConicSolver.get_coeff_offset(con@expr)
      ineq_coeffs[[1]] <- c(ineq_coeffs[[1]], coeff_offset[[1]])
      ineq_coeffs[[2]] <- c(ineq_coeffs[[2]], coeff_offset[[2]])
    }
    Fmat <- Matrix(do.call(rbind, ineq_coeffs[[1]]), byrow = TRUE, sparse = TRUE)
    g <- -ineq_coeffs[[2]]
  } else {
    Fmat <- Matrix(nrow = 0, ncol = n, byrow = TRUE, sparse = TRUE)
    g <- -matrix(nrow = 0, ncol = 0)
  }
  
  # Create dictionary with problem data.
  variables <- variables(problem)[[1]]
  data <- list()
  data[P_KEY] <- Matrix(P, sparse = TRUE)
  data[Q_KEY] <- q
  data[A_KEY] <- Matrix(A, sparse = TRUE)
  data[B_KEY] <- b
  data[F_KEY] <- Matrix(Fmat, sparse = TRUE)
  data[G_KEY] <- g
  data[BOOL_IDX] <- sapply(variables@boolean_idx, function(t) { t[[1]] })
  data[INT_IDX] <- sapply(variables@integer_idx, function(t) { t[[1]] })
  data$n_var <- n
  data$n_eq <- nrow(A)
  data$n_ineq <- nrow(Fmat)
  
  inverse_data@sorted_constraints <- c(eq_cons, ineq_cons)
  
  # Add information about integer variables.
  inverse_data@is_mip <- length(data[BOOL_IDX]) > 0 || length(data[INT_IDX]) > 0
  
  return(list(data, inverse_data))
})
