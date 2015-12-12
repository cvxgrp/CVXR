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
  else if(length(constr_map[[SDP_MAP]]) > 0)
    return(CVXOPT())
  # Otherwise use ECOS
  else
    return(ECOS())
}

setMethod("validate_solver", "Solver", function(solver, constraints) {
  constr_map <- SymData.filter_constraints(constraints)
  if(( (length(constr_map[[BOOL_MAP]]) > 0 || length(constr_map[[INT_MAP]]) > 0) && !mip_capable(solver)) ||
     (length(constr_map[[SDP_MAP]]) > 0 && !sdp_capable(solver)) ||
     (length(constr_map[[EXP_MAP]]) > 0 && !exp_capable(solver)) ||
     (length(constr_map[[SOC_MAP]]) > 0 && !socp_capable(solver)) || 
     (length(constraints) == 0 && name(solver) %in% c(SCS_NAME, GLPK_NAME)))
    stop("The solver ", name(solver), " cannot solve the problem")
})

setMethod("validate_cache", "Solver", function(solver, objective, constraints, cached_data) {
  prob_data <- cached_data[[name(solver)]]
  if(!is.null(prob_data@sym_data) && (objective != prob_data@sym_data@objective || constraints != prob_data@sym_data@constraints)) {
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
  prob_data@sym_data
})

setMethod("get_matrix_data", "Solver", function(solver, objective, constraints, cached_data) {
  sym_data <- get_sym_data(objective, constraints, cached_data)
  prob_data <- cached_data[[name(solver)]]
  if(is.na(prob_data@matrix_data))
    prob_data@matrix_data <- MatrixData(sym_data, matrix_intf(solver), vec_intf(solver), solver)
  prob_data@matrix_data
})

setMethod("get_problem_data", "Solver", function(solver, objective, constraints, cached_data) {
  sym_data <- get_sym_data(objective, constraints, cached_data)
  matrix_data <- get_matrix_data(objective, constraints, cached_data)
  
  data <- list()
  obj <- get_objective(matrix_data)
  eq <- get_eq_constr(matrix_data)
  ineq <- get_ineq_constr(matrix_data)
  
  data[[C]] <- obj[[1]]
  data[[OFFSET]] <- obj[[2]]
  data[[A]] <- eq[[1]]
  data[[B]] <- eq[[2]]
  data[[G]] <- ineq[[1]]
  data[[H]] <- ineq[[2]]
  data[[DIMS]] <- sym_data@dims
  
  conv_idx <- Solver._noncvx_id_to_idx(data[[DIMS]], sym_data@var_offsets, sym_data@var_sizes)
  data[[BOOL_IDX]] <- conv_idx$bool_idx
  data[[INT_IDX]] <- conv_idx$int_idx
  data
})

Solver.is_mip <- function(data) {
  length(data[BOOL_IDX]) > 0 || length(data[INT_IDX]) > 0
}

Solver._noncvx_id_to_idx <- function(dims, var_offsets, var_sizes) {
  bool_idx <- lapply(dims[BOOL_IDS], function(var_id) {
    offset <- var_offsets[var_id]
    size <- var_sizes[var_id]
    offset + seq(1, size[1]*size[2], by = 1)
  })
  
  int_idx <- lapply(dims[INT_IDS], function(var_id) {
    offset <- var_offsets[var_id]
    size <- var_sizes[var_id]
    offset + seq(1, size[1]*size[2], by = 1)
  })
  
  list(bool_idx = bool_idx, int_idx = int_idx)
}

setClass("ECOS", contains = "Solver")
ECOS <- function() { new("ECOS") }

setMethod("name", "ECOS", function(object) { ECOS_NAME })
setMethod("matrix_intf", "ECOS", function(solver) { DEFAULT_SPARSE_INTF })
setMethod("vec_intf", "ECOS", function(solver) { DEFAULT_INTF })

setMethod("split_constr", "ECOS", function(solver, constr_map) {
  list(eq_constr = constr_map[[EQ_MAP]], ineq_constr = constr_map[[LEQ_MAP]], nonlin_constr = list())
})

setMethod("ecos_solve", "ECOS", function(solver, objective, constraints, cached_data, warm_start, verbose, solver_opts) {
  require(ECOSolveR)
  data <- get_problem_data(solver, objective, constraints, cached_data)
  data[[DIMS]]['e'] <- data[[DIMS]][[EXP_DIM]]
  results_dict <- ECOSolveR::ECOS_csolve(c = data[[C]], G = data[[G]], h = data[[H]], dims = data[[DIMS]], A = data[[A]], b = data[[B]])
  format_results(solver, result_dict, data, cached_data)
})

setMethod("format_results", "ECOS", function(solver, results_dict, data, cached_data) {
  new_results <- list()
  status <- status_map(solver, results_dict['info']['exitFlag'])
  new_results[[STATUS]] <- status
  if(new_results[[STATUS]] %in% SOLUTION_PRESENT) {
    primal_val <- results_dict['info']['pcost']
    new_results[[VALUE]] <- primal_val + data[[OFFSET]]
    new_results[[PRIMAL]] <- results_dict['x']
    new_results[[EQ_DUAL]] <- results_dict['y']
    new_results[[INEQ_DUAL]] <- results_dict['z']
  }
  new_results
})
