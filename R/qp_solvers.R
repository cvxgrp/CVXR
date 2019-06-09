# QPSolver requires objectives to be stuffed in the following way.
is_stuffed_qp_objective <- function(objective) {
  expr <- objective@expr
  return(class(expr) == "AddExpression" && length(expr@args) == 2 && class(expr@args[[1]]) == "QuadForm" && class(expr@args[[2]]) == "MulExpression" && is_affine(expr@args[[2]]))
}

# A QP solver interface.
setClass("QpSolver", contains = "ReductionSolver")

setMethod("accepts", signature(object = "QpSolver", problem = "Problem"), function(object, problem) {
  return(class(problem@objective) == "Minimize" && is_stuffed_qp_objective(problem@objective) && are_args_affine(problem@constraints) &&
           all(sapply(problem@constraints, function(c) { class(c) == "ZeroConstraint" || class(c) == "NonPosConstraint" })))
})

setMethod("perform", signature(object = "QpSolver", problem = "Problem"), function(object, problem) {
  # Construct QP problem data stored in a dictionary.
  # The QP has the following form
  #    minimize 1/2 x' P x + q' x
  #    subject to A x = b
  #               F x <= g
  inverse_data <- InverseData(problem)

  obj <- problem@objective
  # quadratic part of objective is t(x) %*% P %*% x, but solvers expect 0.5*t(x) %*% P %*% x.
  P <- 2*value(obj@expr@args[[1]]@args[[2]])
  q <- as.vector(value(obj@expr@args[[2]]@args[[1]]))

  # Get number of variables.
  n <- problem@.size_metrics@num_scalar_variables

  if(length(problem@constraints) == 0) {
    eq_cons <- list()
    ineq_cons <- list()
  } else {
    eq_cons <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "ZeroConstraint" })]
    ineq_cons <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "NonPosConstraint" })]
  }

  # TODO: This dependence on ConicSolver is hacky; something should change here.
  if(length(eq_cons) > 0) {
    eq_coeffs <- list(list(), c())
    for(con in eq_cons) {
      coeff_offset <- ConicSolver.get_coeff_offset(con@expr)
      eq_coeffs[[1]] <- c(eq_coeffs[[1]], list(coeff_offset[[1]]))
      eq_coeffs[[2]] <- c(eq_coeffs[[2]], coeff_offset[[2]])
    }
    A <- Matrix(do.call(rbind, eq_coeffs[[1]]), sparse = TRUE)
    b <- -eq_coeffs[[2]]
  } else {
    A <- Matrix(nrow = 0, ncol = n, sparse = TRUE)
    b <- -matrix(nrow = 0, ncol = 0)
  }

  if(length(ineq_cons) > 0) {
    ineq_coeffs <- list(list(), c())
    for(con in ineq_cons) {
      coeff_offset <- ConicSolver.get_coeff_offset(con@expr)
      ineq_coeffs[[1]] <- c(ineq_coeffs[[1]], list(coeff_offset[[1]]))
      ineq_coeffs[[2]] <- c(ineq_coeffs[[2]], coeff_offset[[2]])
    }
    Fmat <- Matrix(do.call(rbind, ineq_coeffs[[1]]), sparse = TRUE)
    g <- -ineq_coeffs[[2]]
  } else {
    Fmat <- Matrix(nrow = 0, ncol = n, sparse = TRUE)
    g <- -matrix(nrow = 0, ncol = 0)
  }

  # Create dictionary with problem data.
  variables <- variables(problem)[[1]]
  data <- list()
  data[[P_KEY]] <- Matrix(P, sparse = TRUE)
  data[[Q_KEY]] <- q
  data[[A_KEY]] <- Matrix(A, sparse = TRUE)
  data[[B_KEY]] <- b
  data[[F_KEY]] <- Matrix(Fmat, sparse = TRUE)
  data[[G_KEY]] <- g
  data[[BOOL_IDX]] <- sapply(variables@boolean_idx, function(t) { t[[1]] })
  data[[INT_IDX]] <- sapply(variables@integer_idx, function(t) { t[[1]] })
  data$n_var <- n
  data$n_eq <- nrow(A)
  data$n_ineq <- nrow(Fmat)

  inverse_data@sorted_constraints <- c(eq_cons, ineq_cons)

  # Add information about integer variables.
  inverse_data@is_mip <- length(data[[BOOL_IDX]]) > 0 || length(data[[INT_IDX]]) > 0

  return(list(data, inverse_data))
})

# QP interface for the CPLEX solver.
CPLEX_QP <- setClass("CPLEX_QP", contains = "QpSolver")

setMethod("mip_capable", "CPLEX_QP", function(solver) { TRUE })

# Map of CPLEX status to CVXR status.
# TODO: Add more!
setMethod("status_map", "CPLEX_QP", function(solver, status) {
  if(status %in% c(1, 101))
    OPTIMAL
  else if(status %in% c(3, 22, 4, 103))
    INFEASIBLE
  else if(status %in% c(2, 21, 118))
    UNBOUNDED
  else if(status %in% c(10, 107))
    USER_LIMIT
  else
    stop("CPLEX status unrecognized: ", status)
})

setMethod("name", "CPLEX_QP", function(x) { CPLEX_NAME })
setMethod("import_solver", "CPLEX_QP", function(solver) { requireNamespace("Rcplex", quietly = TRUE) })
setMethod("alt_invert", "CPLEX_QP", function(object, results, inverse_data) {
  model <- results$model
  attr <- list()
  if("cputime" %in% names(results))
    attr[SOLVE_TIME] <- results$cputime
  attr[NUM_ITERS] <- as.integer(get_num_barrier_iterations(model@solution@progress))

  status <- status_map(object, get_status(model@solution))

  if(status %in% SOLUTION_PRESENT) {
    # Get objective value.
    opt_val <- get_objective_value(model@solution)

    # Get solution.
    x <- as.matrix(get_values(model@solution))
    primal_vars <- list()
    primal_vars[names(inverse_data@id_map)[1]] <- x

    # Only add duals if not a MIP.
    dual_vars <- NA
    if(!inverse_data$is_mip) {
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

setMethod("solve_via_data", "CPLEX_QP", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
  requireNamespace("Rcplex", quietly = TRUE)
  P <- Matrix(data[[P_KEY]], byrow = TRUE, sparse = TRUE)
  q <- data[[Q_KEY]]
  A <- Matrix(data[[A_KEY]], byrow = TRUE, sparse = TRUE)
  b <- data[[B_KEY]]
  Fmat <- Matrix(data[[F_KEY]], byrow = TRUE, sparse = TRUE)
  g <- data[[G_KEY]]
  n_var <- data$n_var
  n_eq <- data$n_eq
  n_ineq <- data$n_ineq

  # Define CPLEX problem.
  model <- Rcplex::Cplex()

  # Minimize problem.
  model@objective@set_sense(model@objective@sense@minimize)

  # Add variables and linear objective.
  var_idx <- cplex_add(model@variables, obj = q, lb = rep(-Inf, n_var), ub = rep(Inf, n_var))

  # Constraint binary/integer variables if present.
  for(i in data[BOOL_IDX])
    set_types(model@variables, var_idx[i], model@variables@type@binary)

  for(i in data[INT_IDX])
    set_types(model@variables, var_idx[i], model@variables@type@integer)

  # Add constraints.
  lin_expr <- list()
  rhs <- list()
  for(i in 1:n_eq) {   # Add equalities.
    start <- A@p[i]
    end <- A@p[i+1]
    row <- list(A@i[start:end], A@j[start:end], A@x[start:end])
    lin_expr <- c(lin_expr, list(row))
    rhs <- c(rhs, b[i])
  }
  if(length(lin_expr) > 0)
    cplex_add(model@linear_constraints, lin_expr = lin_expr, senses = rep("E", length(lin_expr)), rhs = rhs)

  lin_expr <- list()
  rhs <- list()
  for(i in 1:n_ineq) {   # Add inequalities.
    start <- Fmat@p[i]
    end <- Fmat@p[i+1]
    row <- list(Fmat@i[start:end], Fmat@j[start:end], Fmat@x[start:end])
    lin_expr <- c(lin_expr, list(row))
    rhs <- c(rhs, g[i])
  }
  if(length(lin_expr) > 0)
    cplex_add(model@linear_constraints, lin_expr = lin_expr, senses = rep("L", length(lin_expr)), rhs = rhs)

  # Set quadratic cost.
  if(nnzero(P) > 0) {   # Only if quadratic form is not null.
    qmat <- list()
    for(i in 1:n_var) {
      start <- P@p[i]
      end <- P@p[i+1]
      qmat <- c(qmat, list(P@i[start:end], P@j[start:end], P@x[start:end]))
    }
    set_quadratic(model@objective, qmat)
  }

  # Set parameters.
  if(!verbose) {
    set_results_stream(model, NA)
    set_log_stream(model, NA)
    set_error_stream(model, NA)
    set_warning_stream(model, NA)
  }

  # TODO: The code in CVXR/problems/solvers.R sets CPLEX parameters in the same way,
  # and the code is duplicated here. This should be refactored.
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

  if(!is.null(kwargs) || !is.na(kwargs) || length(kwargs) > 0)
    stop("Invalid keyword argument ", kwargs[[1]])

  # Solve problem.
  results_dict <- list()
  tryCatch({
      start <- get_time(model)
      solve(model)
      end <- get_time(model)
      results_dict$cputime <- end - start
    }, error = function(e) {
      results_dict$status <- SOLVER_ERROR
    }
  )
  results_dict$model <- model
  return(results_dict)
})

# QP interface for the GUROBI solver.
GUROBI_QP <- setClass("GUROBI_QP", contains = "QpSolver")

setMethod("mip_capable", "GUROBI_QP", function(solver) { TRUE })

# Map of GUROBI status to CVXR status.
setMethod("status_map", "GUROBI_QP", function(solver, status) {
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

setMethod("name", "GUROBI_QP", function(x) { GUROBI_NAME })
setMethod("import_solver", "GUROBI_QP", function(solver) {
  requireNamespace("gurobi", quietly = TRUE)
})

##DK: IS THIS FUNCTION NECESSARY ANYMORE WITH invert?
setMethod("alt_invert", "GUROBI_QP", function(object, results, inverse_data) {
  model <- results$model
  x_grb <- getVars(model)
  n <- length(x_grb)
  constraints_grb <- getConstrs(model)
  m <- length(constraints_grb)

  # Start populating attribute dictionary.
  attr <- list()
  attr[SOLVE_TIME] <- model@Runtime
  attr[NUM_ITERS] <- model@BarIterCount

  # Map GUROBI statuses back to CVXR statuses.
  status <- status_map(object, model@Status)

  if(status %in% SOLUTION_PRESENT) {
    opt_val <- model@objVal
    x <- as.matrix(sapply(1:n, function(i) { x_grb[i]@X }))

    primal_vars <- list()
    primal_vars[names(inverse_data@id_map)[1]] <- x

    # Only add duals if not a MIP.
    dual_vars <- NA
    if(!inverse_data$is_mip) {
      y <- -as.matrix(sapply(1:m, function(i) { constraints_grb[i]@Pi }))
      dual_vars <- get_dual_values(y, extract_dual_value, inverse_data@sorted_constraints)
    } else {
      primal_vars <- NA
      dual_vars <- NA
      opt_val <- Inf
      if(status == UNBOUNDED)
        opt_val <- -Inf
    }
  }

  return(Solution(status, opt_val, primal_vars, dual_vars, attr))
})

setMethod("solve_via_data", "GUROBI_QP", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
  requireNamespace("gurobi", quietly = TRUE)
  # N.B. Here we assume that the matrices in data are in CSC format.
  P <- data[[P_KEY]]   # TODO: Convert P matrix to COO format?
  q <- data[[Q_KEY]]
  A <- data[[A_KEY]]   # TODO: Convert A matrix to CSR format?
  b <- data[[B_KEY]]
  Fmat <- data[[F_KEY]]   # TODO: Convert F matrix to CSR format?
  g <- data[[G_KEY]]
  n <- data$n_var

  # Create a new model.
  model <- list()

  #Doesn't seem like adding variables exists in the R version, but setting var type is
  # Add variables.
  
  #Add variable types
  vtype <- character(n)
  for(i in seq_along(data[[BOOL_IDX]])){
    vtype[data[[BOOL_IDX]][i]] <- 'B' #B for binary
  }
  for(i in seq_along(data[[INT_IDX]])){
    vtype[data[[INT_IDX]][i]] <- 'I' #I for integer
  }
  
  for(i in 1:n) {
    if(vtype[i] == ""){
      vtype[i] <- 'C' #C for continuous
    }
    #addVar(model, ub = Inf, lb = -Inf, vtype = vtype): don't need in R
  }
  model$vtype <- vtype #put in variable types
  model$lb <- rep(-Inf, n)
  model$ub <- rep(Inf, n)
  
  #update(model): doesn't exist in R
  #x <- getVars(model)

  # Add equality constraints: iterate over the rows of A,
  # adding each row into the model.
  model$A <- rbind(A, Fmat)
  model$rhs <- c(b, g)
  model$sense <- c(rep('=', dim(A)[1]), rep('<', dim(Fmat)[1])) #in their documentation they just have < than not <=?
  
  # if(!is.null(model$v811_addMConstrs)) {
  #   # @e can pass all of A == b at once.
  #   sense <- rep(GRB.EQUAL, nrow(A))
  #   model <- model$v811_addMConstrs(A, sense, b)              What is the R equivalent of this??
  # } else if(nrow(A) > 0) {
  #   for(i in 1:nrow(A)) {                                     Is this bit ever necessary in R?
  #     start <- A@p[i]
  #     end <- A@p[i+1]
  #     # variables <- mapply(function(i, j) { x[i,j] }, A@i[start:end], A@j[start:end])   # Get nnz.
  #     variables <- x[A@i[start:end]]
  #     coeff <- A@x[start:end]
  #     expr <- gurobi::LinExpr(coeff, variables)
  #     addConstr(model, expr, "equal", b[i])
  #   }
  # }
  # update(model)

  # Add inequality constraints: iterate over the rows of F,
  # adding each row into the model.
  # if(!is.null(model$v811_addMConstrs)) {
  #   # We can pass all of F <= g at once.
  #   sense <- rep(GRB.LESS_EQUAL, nrow(Fmat))
  #   model <- model$v811_addMConstrs(Fmat, sense, g)
  # } else if(nrow(Fmat) > 0) {
  #   for(i in nrow(Fmat)) {
  #     start <- Fmat@p[i]
  #     end <- Fmat@p[i+1]
  #     # variables <- mapply(function(i, j) { x[i,j] }, Fmat@i[start:end], Fmat@j[start:end])   # Get nnz.
  #     variables <- x[Fmat@i[start:end]]
  #     coeff <- Fmat@x[start:end]
  #     expr <- gurobi::LinExpr(coeff, variables)
  #     addConstr(model, expr, "less_equal", g[i])
  #   }
  # }
  #update(model) Not in R

  # Define objective.
  
  
  ####CHECK MATH####
  #Conjecture P is Q matrix and q is c, which is obj for gurobi
  model$Q <- P*.5
  model$obj <- q 
  
  # obj <- gurobi::QuadExpr()
  # if(!is.null(model$v811_setMObjective))
  #   model <- model$v811_setMObjective(0.5*P, q)
  # else {
  #   nnz <- nnzero(P)
  #   if(nnz > 0) {   # If there are any nonzero elements in P.
  #     for(i in 1:nnz)
  #       gurobi_add(obj, 0.5*P@x[i]*x[P@i[i]]*x[P@j[i]])
  #   }
  #   gurobi_add(obj, gurobi::LinExpr(q, x))   # Add linear part.
  #   setObjective(model, obj)   # Set objective.
  # }
  # update(model)

  # Set verbosity and other parameters.
  params <- list()
  #setParam(model, "OutputFlag", verbose)
  params$OutputFlag <- as.numeric(verbose)
  # TODO: User option to not compute duals.
  params$QCPDual <- 1 #equivalent to TRUE

  # for(key in names(solver_opts))
  #   setParam(model, key, solver_opts[key])
  for(i in seq_along(solver_opts)){
    params[[ names(solver_opts)[i] ]] <- solver_opts[i] 
  }
  
  # Update model. Not a thing in R
  #update(model)

  # Solve problem.
  #results_dict <- gurobi(model, params)
  results_dict <- list()
  tryCatch({
    results_dict$solution <- gurobi::gurobi(model, params)   # Solve.
  }, error = function(e) {   # Error in the solution.
    results_dict$status <- 'SOLVER_ERROR'
  })
  results_dict$model <- model
  
  return(results_dict)
})

#DK WRITTEN FUNCTION
setMethod("invert", signature(object = "GUROBI_QP", solution = "list", inverse_data = "InverseData"), function(object, solution, inverse_data){
  model <- solution$model
  solution <- solution$solution
  x_grb <- model$x
  n <- length(x_grb)
  constraints_grb <- model$rhs
  m = length(constraints_grb)
  
  attr <- list()
  attr[[SOLVE_TIME]] <- solution$runtime
  attr[[NUM_ITERS]] <- solution$baritercount
  
  status <- status_map(object, solution$status)
  
  if(status %in% SOLUTION_PRESENT){
    opt_val <- solution$objval
    x <- solution$x
    
    primal_vars <- list()
    primal_vars[[names(inverse_data@id_map)[1]]] <- x
    
    #Only add duals if not a MIP
    dual_vars <- NA
    if(!inverse_data@is_mip){
      y <- solution$pi
      dual_vars <- get_dual_values(y, extract_dual_value, inverse_data@sorted_constraints)
    } else{
      primal_vars <- NA
      dual_vars <- NA
      opt_val <- Inf
      if(status == UNBOUNDED){
        opt_val <- -Inf
      }
    }
  }
  return(Solution(status, opt_val, primal_vars, dual_vars, attr))
})

# QP interface for the OSQP solver.
OSQP <- setClass("OSQP", contains = "QpSolver")

# Map of OSQP status to CVXPY status.
setMethod("status_map", "OSQP", function(solver, status) {
  if(status == 1)
    OPTIMAL
  else if(status == 2)
    OPTIMAL_INACCURATE
  else if(status == -3)
    INFEASIBLE
  else if(status == 3)
    INFEASIBLE_INACCURATE
  else if(status == -4)
    UNBOUNDED
  else if(status == 4)
    UNBOUNDED_INACCURATE
  else if(status == -2 || status == -5 || status == -10)   # -2: Maxiter reached. -5: Interrupted by user. -10: Unsolved.
    SOLVER_ERROR
 else
    stop("OSQP status unrecognized: ", status)
})

setMethod("name", "OSQP", function(x) { OSQP_NAME })
setMethod("import_solver", "OSQP", function(solver) {
  requireNamespace("osqp", quietly = TRUE)
})

## DWK CHANGE
## solution = "Solution" before. CVXPY has osqp.OSQP.results class
## In R, OSQP returns a list, so I changed the solution signature into a list
## setMethod("invert", signature(object = "OSQP", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
## DWK CHANGE END
setMethod("invert", signature(object = "OSQP", solution = "list", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  attr <- list()
  ## DWK CHANGES Below
  ## Change slots to list elements
  attr[[SOLVE_TIME]] <- solution$info$run_time

  # Map OSQP statuses back to CVXR statuses.
  status <- status_map(object, solution$info$status_val)
  ## DWK CHANGE END
  if(status %in% SOLUTION_PRESENT) {
      ## DWK CHANGES Below
      opt_val <- solution$info$obj_val
      ## DWK CHANGE END
    primal_vars <- list()
    ## DWK CHANGE
    ## primal_vars[names(inverse_data@id_map)] <- solution@x
    primal_vars[[names(inverse_data@id_map)]] <- solution$x
    dual_vars <- get_dual_values(solution$y, extract_dual_value, inverse_data@sorted_constraints)
    attr[[NUM_ITERS]] <- solution$info$iter
    return(Solution(status, opt_val, primal_vars, dual_vars, attr))
    ## DWK CHANGE END
  } else
    return(failure_solution(status))
})

setMethod("solve_via_data", "OSQP", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
  requireNamespace("osqp", quietly = TRUE)
  P <- data[[P_KEY]]
  q <- data[[Q_KEY]]
  A <- Matrix(do.call(rbind, list(data[[A_KEY]], data[[F_KEY]])), sparse = TRUE)
  data$full_A <- A
  uA <- c(data[[B_KEY]], data[[G_KEY]])
  data$u <- uA
  lA <- c(data[[B_KEY]], rep(-Inf, length(data[[G_KEY]])))
  data$l <- lA

  # Overwrite defaults eps_abs = eps_rel = 1e-3, max_iter = 4000.
  if(is.null(solver_opts$eps_abs))
    solver_opts$eps_abs <- 1e-4
  if(is.null(solver_opts$eps_rel))
    solver_opts$eps_rel <- 1e-4
  if(is.null(solver_opts$max_iter))
    solver_opts$max_iter <- 10000

  if(!is.null(solver_cache) && length(solver_cache) > 0 && name(object) %in% names(solver_cache)) {
    # Use cached data.
    cache <- solver_cache[[name(object)]]
    solver <- cache[[1]]
    old_data <- cache[[2]]
    results <- cache[[3]]

    same_pattern <- all(dim(P) == dim(old_data[[P_KEY]])) && all(P@p == old_data[[P_KEY]]@p) && all(P@i == old_data[[P_KEY]]@i) &&
                    all(dim(A) == dim(old_data$full_A)) && all(A@p == old_data$full_A@p) && all(A@i == old_data$full_A@i)
  } else
    same_pattern <- FALSE

  # TODO: Check syntax for setting up R's OSQP solver.
  # If sparsity pattern differs, need to do setup.
  if(warm_start && same_pattern) {
    new_args <- list()
    for(key in c("q", "l", "u")) {
      if(any(data[[key]] != old_data[[key]]))
        new_args[[key]] <- data[[key]]
    }
    factorizing <- FALSE

    if(any(P@x != old_data[[P_KEY]]@x)) {
      P_triu <- triu(P)
      new_args$Px <- P_triu@x
      factorizing <- TRUE
    }
    if(any(A@x != old_data$full_A@x)) {
      new_args$Ax <- A@x
      factorizing <- TRUE
    }

    if(length(new_args) > 0)
      update(solver, new_args)

    # Map OSQP statuses back to CVXR statuses.
    status <- status_map(object, results@info@status_val)
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
      ## DWK CHANGE
      ## solver <- osqp::OSQP()
      solver <- osqp::osqp(P, q, A, lA, uA, solver_opts)
      ## setup(solver, P, q, A, lA, uA, verbose = verbose, solver_opts)
      ## DWK CHANGE END


  }

  ## DWK CHANGE
  ## results <- solve(solver)
  results <- solver$Solve()
  ## if(!is.null(solver_cache))
  if(identical(solver_cache, list()) || !is.null(solver_cache))   # solver_cache is a list() object, so it throws an error here
      solver_cache[[name(object)]] <- list(solver, data, results)
  ## DWK CHANGE END
  return(results)
})
