Solver <- setClass("Solver", contains = "VIRTUAL")

Solver.choose_solver <- function(constraints) {
  constr_map <- SymData.filter_constraints(constraints)
  # If no constraints, use ECOS.
  if(length(constraints) == 0)
    return(ECOS())
  # If mixed integer constraints, use ECOS_BB.
  else if(length(constr_map[[BOOL_MAP]]) > 0 || length(constr_map[[INT_MAP]]) > 0)
    return(ECOS_BB())
  # If SDP, defaults to CVXOPT.
  else if(length(constr_map[[SDP_MAP]]) > 0) {
    if("cvxopt" %in% installed.packages()) {
      require(cvxopt)
      return(CVXOPT())
    } else
      return(SCS())
  }
  # Otherwise use ECOS
  else
    return(ECOS())
}

setMethod("validate_solver", "Solver", function(solver, constraints) {
  # Check the solver is installed.
  if(!import_solver(solver))
    stop("The solver ", name(solver), " is not installed.")

  # Check the solver can solve the problem.
  constr_map <- SymData.filter_constraints(constraints)

  if( (length(constr_map[[BOOL_MAP]]) > 0 || length(constr_map[[INT_MAP]]) > 0) && !mip_capable(solver))
    Solver._reject_problem(solver, "it cannot solve mixed-integer problems")
  else if(length(constr_map[[SDP_MAP]]) > 0 && !sdp_capable(solver))
    Solver._reject_problem(solver, "it cannot solve semidefinite problems")
  else if(length(constr_map[[EXP_MAP]]) > 0 && !exp_capable(solver))
    Solver._reject_problem(solver, "it cannot solve exponential cone problems")
  else if(length(constr_map[[SOC_MAP]]) > 0 && !socp_capable(solver))
    Solver._reject_problem(solver, "it cannot solve second-order cone problems")
  else if(length(constraints) == 0 && name(solver) %in% c("SCS", "GLPK"))
    Solver._reject_problem(solver, "it cannot solve unconstrained problems")
})

Solver._reject_problem <- function(solver, reason) {
  message <- paste("The solver", name(solver), "cannot solve the problem because", reason)
  stop(message)
}

setMethod("validate_cache", "Solver", function(solver, objective, constraints, cached_data) {
  prob_data <- cached_data[[name(solver)]]
  if(!is.null(prob_data@sym_data) && (!isTRUE(all.equal(objective, prob_data@sym_data@objective)) ||
                                      !isTRUE(all.equal(constraints, prob_data@sym_data@constraints)))) {
    prob_data@sym_data <- NULL
    prob_data@matrix_data <- NULL
  }
  cached_data[[name(solver)]] <- prob_data
  cached_data
})

setMethod("get_sym_data", "Solver", function(solver, objective, constraints, cached_data) {
  cached_data <- validate_cache(solver, objective, constraints, cached_data)
  prob_data <- cached_data[[name(solver)]]
  if(is.null(prob_data@sym_data))
    prob_data@sym_data <- SymData(objective, constraints, solver)
  cached_data[[name(solver)]] <- prob_data
  cached_data
})

setMethod("get_matrix_data", "Solver", function(solver, objective, constraints, cached_data) {
  cached_data <- get_sym_data(solver, objective, constraints, cached_data)
  sym_data <- cached_data[[name(solver)]]@sym_data
  prob_data <- cached_data[[name(solver)]]
  if(is.null(prob_data@matrix_data))
    prob_data@matrix_data <- MatrixData(sym_data, solver)   # TODO: Update this constructor for nonlinear constraints
  cached_data[[name(solver)]] <- prob_data
  cached_data
})

setMethod("Solver.get_problem_data", "Solver", function(solver, objective, constraints, cached_data) {
  cached_data <- get_sym_data(solver, objective, constraints, cached_data)
  sym_data <- cached_data[[name(solver)]]@sym_data
  cached_data <- get_matrix_data(solver, objective, constraints, cached_data)
  matrix_data <- cached_data[[name(solver)]]@matrix_data

  data <- list()
  obj <- get_objective(matrix_data)
  eq <- get_eq_constr(matrix_data)
  ineq <- get_ineq_constr(matrix_data)

  data[[C_KEY]] <- obj[[1]]
  data[[OFFSET]] <- obj[[2]]
  data[[A_KEY]] <- eq[[1]]
  data[[B_KEY]] <- eq[[2]]
  data[[G_KEY]] <- ineq[[1]]
  data[[H_KEY]] <- ineq[[2]]
  data[[DIMS]] <- sym_data@.dims

  conv_idx <- Solver._noncvx_id_to_idx(data[[DIMS]], sym_data@.var_offsets, sym_data@.var_sizes)
  data[[DIMS]] <- conv_idx$dims
  data[[BOOL_IDX]] <- conv_idx$bool_idx
  data[[INT_IDX]] <- conv_idx$int_idx
  data
})

setMethod("nonlin_constr", "Solver", function(solver) { FALSE })

Solver.is_mip <- function(data) {
  length(data[BOOL_IDX]) > 0 || length(data[INT_IDX]) > 0
}

Solver._noncvx_id_to_idx <- function(dims, var_offsets, var_sizes) {
  if(BOOL_IDS %in% names(dims)) {
    bool_idx <- lapply(dims[[BOOL_IDS]], function(var_id) {
      offset <- var_offsets[var_id]
      size <- var_sizes[var_id]
      offset + seq(1, size[1]*size[2], by = 1)
    })
    dims[[BOOL_IDS]] <- NULL
  } else
    bool_idx <- list()

  if(INT_IDS %in% names(dims)) {
    int_idx <- lapply(dims[[INT_IDS]], function(var_id) {
      offset <- var_offsets[var_id]
      size <- var_sizes[var_id]
      offset + seq(1, size[1]*size[2], by = 1)
    })
    dims[[INT_IDS]] <- NULL
  } else
    int_idx <- list()

  list(dims = dims, bool_idx = bool_idx, int_idx = int_idx)
}

setClass("ECOS", contains = "Solver")
ECOS <- function() {
    new("ECOS")
    ##ECOS$new()
}

# ECOS capabilities
setMethod("lp_capable", "ECOS", function(solver) { TRUE })
setMethod("socp_capable", "ECOS", function(solver) { TRUE })
setMethod("sdp_capable", "ECOS", function(solver) { FALSE })
setMethod("exp_capable", "ECOS", function(solver) { TRUE })
setMethod("mip_capable", "ECOS", function(solver) { FALSE })

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

# Map of ECOS status to CVXPY status.
setMethod("status_map", "ECOS", function(solver, status) {
  if(status == 0) OPTIMAL
  else if(status == 1) INFEASIBLE
  else if(status == 2) UNBOUNDED
  else if(status == 10) OPTIMAL_INACCURATE
  else if(status == 11) INFEASIBLE_INACCURATE
  else if(status == 12) UNBOUNDED_INACCURATE
  else if(status %in% c(-1, -2, -3, -4, -7)) SOLVER_ERROR
  else stop("ECOS status unrecognized: ", status)
})

setMethod("name", "ECOS", function(object) { ECOS_NAME })
setMethod("import_solver", "ECOS", function(solver) { require(ECOSolveR) })
setMethod("matrix_intf", "ECOS", function(solver) { DEFAULT_SPARSE_INTF })
setMethod("vec_intf", "ECOS", function(solver) { DEFAULT_INTF })
setMethod("split_constr", "ECOS", function(solver, constr_map) {
  list(eq_constr = constr_map[[EQ_MAP]], ineq_constr = constr_map[[LEQ_MAP]], nonlin_constr = list())
})

setMethod("Solver.solve", "ECOS", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {
  require(ECOSolveR)
  data <- Solver.get_problem_data(solver, objective, constraints, cached_data)
  data[[DIMS]]['e'] <- data[[DIMS]][[EXP_DIM]]
  
  # TODO: Naras please fix this by making ECOSolveR handle type conversion, e.g. logical -> integer
  if(prod(dim(data[[G_KEY]])) == 0) data[[G_KEY]] <- NULL
  if(prod(dim(data[[A_KEY]])) == 0) data[[A_KEY]] <- NULL
  data[[DIMS]] <- lapply(data[[DIMS]], function(dim) { as.integer(dim) })
  solver_opts <- ECOSolveR::ecos.control()
  solver_opts$verbose <- as.integer(verbose)
  other_opts <- list(...)
  solver_opts[names(other_opts)] <- other_opts
  
  results_dict <- ECOSolveR::ECOS_csolve(c = data[[C_KEY]], G = data[[G_KEY]], h = data[[H_KEY]], dims = data[[DIMS]], A = data[[A_KEY]], b = data[[B_KEY]], control = solver_opts)
  format_results(solver, results_dict, data, cached_data)
})

setMethod("format_results", "ECOS", function(solver, results_dict, data, cached_data) {
  new_results <- list()
  status <- status_map(solver, results_dict$retcodes[["exitFlag"]])
  new_results[[STATUS]] <- status

  # Timing data
  new_results[[SOLVE_TIME]] <- results_dict$timing[["tsolve"]]
  new_results[[SETUP_TIME]] <- results_dict$timing[["tsetup"]]
  new_results[[NUM_ITERS]] <- results_dict$retcodes[["iter"]]

  if(new_results[[STATUS]] %in% SOLUTION_PRESENT) {
    primal_val <- results_dict$summary[["pcost"]]
    new_results[[VALUE]] <- primal_val + data[[OFFSET]]
    new_results[[PRIMAL]] <- results_dict$x
    new_results[[EQ_DUAL]] <- results_dict$y
    new_results[[INEQ_DUAL]] <- results_dict$z
  }
  new_results
})

setClass("SCS", contains = "ECOS")
SCS <- function() {
    new("SCS")
    ##SCS$new()
}

# SCS capabilities
setMethod("lp_capable", "SCS", function(solver) { TRUE })
setMethod("socp_capable", "SCS", function(solver) { TRUE })
setMethod("sdp_capable", "SCS", function(solver) { TRUE })
setMethod("exp_capable", "SCS", function(solver) { TRUE })
setMethod("mip_capable", "SCS", function(solver) { FALSE })

# Map of SCS status to CVXPY status.
setMethod("status_map", "SCS", function(solver, status) {
  if(status == "Solved") OPTIMAL
  else if(status == "Solved/Inaccurate") OPTIMAL_INACCURATE
  else if(status == "Unbounded") UNBOUNDED
  else if(status == "Unbounded/Inaccurate") UNBOUNDED_INACCURATE
  else if(status == "Infeasible") INFEASIBLE
  else if(status == "Infeasible/Inaccurate") INFEASIBLE_INACCURATE
  else if(status %in% c("Failure", "Indeterminate")) SOLVER_ERROR
  else stop("SCS status unrecognized: ", status)
})

setMethod("name", "SCS", function(object) { SCS_NAME })
setMethod("import_solver", "SCS", function(solver) { require(scs) })
setMethod("split_constr", "SCS", function(solver, constr_map) {
  list(eq_constr = c(constr_map[[EQ_MAP]], constr_map[[LEQ_MAP]]), ineq_constr = list(), nonlin_constr = list())
})

setMethod("Solver.solve", "SCS", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) {
  require(scs)
  data <- Solver.get_problem_data(solver, objective, constraints, cached_data)

  # Set the options to be VERBOSE plus any user-specific options
  solver_opts <- list(...)
  solver_opts$verbose <- verbose
  scs_args <- list(A = data[[A_KEY]], b = data[[B_KEY]], obj = data[[C_KEY]], cone = data[[DIMS]])

  # If warm starting, add old primal and dual variables
  solver_cache <- cached_data[name(solver)]
  if(warm_start && !is.na(solver_cache$prev_result)) {
    stop("Warm start currently unimplemented")
    scs_args$x <- solver_cache$prev_result$x
    scs_args$y <- solver_cache$prev_result$y
    scs_args$s <- solver_cache$prev_result$s
  }

  results_dict <- do.call(scs::scs, c(scs_args, list(control = solver_opts)))
  # if(length(solver_opts) == 0)
  #  results_dict <- scs::scs(A = data[[A_KEY]], b = data[[B_KEY]], obj = data[[C_KEY]], cone = data[[DIMS]])
  # else
  #  results_dict <- scs::scs(A = data[[A_KEY]], b = data[[B_KEY]], obj = data[[C_KEY]], cone = data[[DIMS]], control = solver_opts)
  format_results(solver, results_dict, data, cached_data)
})

setMethod("format_results", "SCS", function(solver, results_dict, data, cached_data) {
  solver_cache <- cached_data[name(solver)]
  dims <- data[[DIMS]]
  new_results <- list()
  status <- status_map(solver, results_dict$info$status)
  new_results[[STATUS]] <- status

  # Timing and iteration data
  new_results[[SOLVE_TIME]] <- results_dict$info$solveTime/1000
  new_results[[SETUP_TIME]] <- results_dict$info$setupTime/1000
  new_results[[NUM_ITERS]] <- results_dict$info$iter

  if(new_results[[STATUS]] %in% SOLUTION_PRESENT) {
    # Save previous result for possible future warm start
    solver_cache$prev_result <- list(x = results_dict$x, y = results_dict$y, s = results_dict$s)
    primal_val <- results_dict$info$pobj
    new_results[[VALUE]] <- primal_val + data[[OFFSET]]
    new_results[[PRIMAL]] <- results_dict$x
    new_results[[EQ_DUAL]] <- results_dict$y[1:dims[[EQ_DIM]]]

    y <- results_dict$y[(dims[[EQ_DIM]]+1):length(results_dict$y)]
    if(is.null(dims[[SDP_DIM]])) {
      old_sdp_sizes <- 0
      new_sdp_sizes <- 0
    } else {
      old_sdp_sizes <- sum(floor(dims[[SDP_DIM]] * (dims[[SDP_DIM]] + 1)/2))
      new_sdp_sizes <- sum(dims[[SDP_DIM]]^2)
    }
    if(is.null(dims[[SOC_DIM]]))
      y_offset <- dims[[LEQ_DIM]]
    else
      y_offset <- dims[[LEQ_DIM]] + sum(dims[[SOC_DIM]])
    y_true <- rep(0, length(y) + (new_sdp_sizes - old_sdp_sizes))
    y_true_offset <- y_offset
    if(y_true_offset > 0)
      y_true[1:y_true_offset] <- y[1:y_offset]

    # Expand SDP duals from lower triangular to full matrix, scaling off diagonal entries by 1/sqrt(2)
    for(n in dims[[SDP_DIM]]) {
      if(n > 0) {
        tri <- y[(y_offset + 1):(y_offset + floor(n*(n+1)/2))]
        y_true[(y_true_offset + 1):(y_true_offset + n^2)] <- SCS.tri_to_full(tri, n)
        y_true_offset <- y_true_offset + n^2
        y_offset <- y_offset + floor(n*(n+1)/2)
      }
    }

    if(length(y_true) > y_true_offset)
      y_true[(y_true_offset + 1):length(y_true)] <- y[(y_offset + 1):length(y)]
    new_results[[INEQ_DUAL]] <- y_true
  } else {
    # No result to save
    solver_cache$prev_result <- NULL
  }
  return(new_results)
})

SCS.tri_to_full <- function(lower_tri, n) {
  # Expands floor(n*(n+1)/2) lower triangular to full matrix, with off-diagonal entries scaled by 1/sqrt(2)
  full <- matrix(0, nrow = n, ncol = n)
  for(col in 1:n) {
    for(row in col:n) {
      idx <- row - col + floor(n*(n+1)/2) - floor((n-col)*(n-col+1)/2)
      if(row != col) {
        full[row, col] <- lower_tri[idx]/sqrt(2)
        full[col, row] <- lower_tri[idx]/sqrt(2)
      } else
        full[row, col] <- lower_tri[idx]
    }
  }
  return(matrix(full, nrow = n^2))
}

setClass("ECOS_BB", contains = "ECOS")
ECOS_BB <- function() {
  stop("Unimplemented solver")
  new("ECOS_BB")
  ##ECOS_BB$new()
}

setClass("CVXOPT", contains = "Solver")
CVXOPT <- function() {
  stop("Unimplemented solver")
  new("CVXOPT")
  ##CVXOPT$new()
}

#'
#' Solver utilities
#'
# solver_intf <- list(ECOS(), ECOS_BB(), CVXOPT(), GLPK(), GLPK_MI(), CBC(), SCS(), GUROBI(), Elemental(), MOSEK(), LS())
solver_intf <- list(ECOS(), SCS())
SOLVERS <- solver_intf
names(SOLVERS) <- sapply(solver_intf, function(solver) { name(solver) })

installed_solvers <- function() {
  installed <- list()
  for(i in 1:length(SOLVERS)) {
    if(is_installed(SOLVERS[i]))
      installed <- c(installed, names(SOLVERS)[i])
  }
  installed
}
