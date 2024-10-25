#' An interface to the CBC solver
#'
#' @name CBC_CONIC-class
#' @aliases CBC_CONIC
#' @rdname CBC_CONIC-class
#' @export
setClass("CBC_CONIC",
         prototype = list(MIP_CAPABLE = TRUE,
                          SUPPORTED_CONSTRAINTS = supported_constraints(ConicSolver()),
                          MI_SUPPORTED_CONSTRAINTS = supported_constraints(ConicSolver())),
         contains = "SCS")

#' @rdname CBC_CONIC-class
#' @export
CBC_CONIC <- function() { new("CBC_CONIC") }

#' @param solver,object,x A \linkS4class{CBC_CONIC} object.
#' @param status A status code returned by the solver.
#' @describeIn CBC_CONIC Converts status returned by the CBC solver to its respective CVXPY status.
setMethod("status_map", "CBC_CONIC", function(solver, status) {
  if(status$is_proven_optimal)
    OPTIMAL
  else if(status$is_proven_dual_infeasible || status$is_proven_infeasible)
    INFEASIBLE
  else
    SOLVER_ERROR # probably need to check this the most
})

# Map of CBC_CONIC MIP/LP status to CVXR status.
#' @describeIn CBC_CONIC Converts status returned by the CBC solver to its respective CVXPY status for mixed integer problems.
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

#' @describeIn CBC_CONIC Converts status returned by the CBC solver to its respective CVXPY status for linear problems.
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

#' @describeIn CBC_CONIC Returns the name of the solver
setMethod("name", "CBC_CONIC", function(x) { CBC_NAME })

#' @describeIn CBC_CONIC Imports the solver
setMethod("import_solver", "CBC_CONIC", function(solver) {
    requireNamespace("rcbc", quietly = TRUE)
})

#' @param problem A \linkS4class{Problem} object.
#' @describeIn CBC_CONIC Can CBC_CONIC solve the problem?
setMethod("accepts", signature(object = "CBC_CONIC", problem = "Problem"), function(object, problem) {
  # Can CBC_CONIC solve the problem?
  # TODO: Check if the matrix is stuffed.
  if(!is_affine(problem@objective@args[[1]]))
    return(FALSE)
  for(constr in problem@constraints) {
    if(!inherits(constr, supported_constraints(object)))
      return(FALSE)
    for(arg in constr@args) {
      if(!is_affine(arg))
        return(FALSE)
    }
  }
  return(TRUE)
})

#' @describeIn CBC_CONIC Returns a new problem and data for inverting the new solution.
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

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A list containing data necessary for the inversion.
#' @describeIn CBC_CONIC Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "CBC_CONIC", solution = "list", inverse_data = "list"),  function(object, solution, inverse_data) {
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

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance.
#' @param reltol The relative tolerance.
#' @param abstol The absolute tolerance.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn CBC_CONIC Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "CBC_CONIC", function(object, data, warm_start, verbose, feastol, reltol, abstol, num_iter, solver_opts, solver_cache) {
  if (missing(solver_cache)) solver_cache <- new.env(parent=emptyenv())
  cvar <- data$c
  b <- data$b
  ## Conversion below forced by changes in Matrix package version 1.3.x
  A <- as(as(data$A, "CsparseMatrix"), "dgTMatrix")
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

  #Warnigs for including parameters
  if(!all(c(is.null(feastol), is.null(reltol), is.null(abstol), is.null(num_iter)))) {
    warning("Ignoring inapplicable parameter feastol/reltol/abstol for CBC.")
  }

  result <- rcbc::cbc_solve(
    obj = cvar,
    mat = A,
    row_ub = row_ub,
    row_lb = row_lb,
    col_lb = var_lb,
    col_ub = var_ub,
    is_integer = is_integer,
    max = FALSE,
    cbc_args = solver_opts
  )

  return(list(result))
})

#' An interface for the CPLEX solver
#'
#' @name CPLEX_CONIC-class
#' @aliases CPLEX_CONIC
#' @rdname CPLEX_CONIC-class
#' @export
CPLEX_CONIC <- setClass("CPLEX_CONIC",
                        prototype = list(MIP_CAPABLE = TRUE,
                                         SUPPORTED_CONSTRAINTS = c(supported_constraints(ConicSolver()), "SOC"),
                                         MI_SUPPORTED_CONSTRAINTS = c(supported_constraints(ConicSolver()), "SOC")),
                        contains = "SCS")

#' @rdname CPLEX_CONIC-class
#' @export
CPLEX_CONIC <- function() { new("CPLEX_CONIC") }

#' @param solver,object,x A \linkS4class{CPLEX_CONIC} object.
#' @describeIn CPLEX_CONIC Returns the name of the solver.
setMethod("name", "CPLEX_CONIC", function(x) { CPLEX_NAME })

#' @describeIn CPLEX_CONIC Imports the solver.
setMethod("import_solver", "CPLEX_CONIC", function(solver) { requireNamespace("Rcplex", quietly = TRUE) })

#' @param problem A \linkS4class{Problem} object.
#' @describeIn CPLEX_CONIC Can CPLEX solve the problem?
setMethod("accepts", signature(object = "CPLEX_CONIC", problem = "Problem"), function(object, problem) {
  # Can CPLEX solve the problem?
  # TODO: Check if the matrix is stuffed.
  if(!is_affine(problem@objective@args[[1]]))
    return(FALSE)
  for(constr in problem@constraints) {
    if(!inherits(constr, supported_constraints(object)))
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
#' @param status A status code returned by the solver.
#' @describeIn CPLEX_CONIC Converts status returned by the CPLEX solver to its respective CVXPY status.
setMethod("status_map", "CPLEX_CONIC", function(solver, status) {
  if(status %in% c(1, 101, 102)){
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

#' @describeIn CPLEX_CONIC Returns a new problem and data for inverting the new solution.
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

  # Include Boolean and Integer indices
  variables <- variables(problem)[[1]]
  data[[BOOL_IDX]] <- lapply(variables@boolean_idx, function(t) { t[1] })
  data[[INT_IDX]] <- lapply(variables@integer_idx, function(t) { t[1] })
  inv_data$is_mip <- length(data[[BOOL_IDX]]) > 0 || length(data[[INT_IDX]]) > 0

  return(list(object, data, inv_data))
})

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A list containing data necessary for the inversion.
#' @describeIn CPLEX_CONIC Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "CPLEX_CONIC", solution = "list", inverse_data = "list"), function(object, solution, inverse_data) {
  # Returns the solution to the original problem given the inverse_data.
  model <- solution$model

  status <- status_map(object, model$status)
  primal_vars <- list()
  dual_vars <- list()

  if(status %in% SOLUTION_PRESENT) {
    #Get objective value
    opt_val <- model$obj + inverse_data[[OFFSET]]

    #Get solution
    primal_vars[[as.character(inverse_data[[object@var_id]])]] <- model$xopt

    #Only add duals if not a MIP
    if(!inverse_data[["is_mip"]]) {
      #eq and leq constraints all returned at once by CPLEX
      eq_dual <- get_dual_values(solution$eq_dual, extract_dual_value, inverse_data[[object@eq_constr]])
      leq_dual <- get_dual_values(solution$ineq_dual, extract_dual_value, inverse_data[[object@neq_constr]])
      eq_dual <- utils::modifyList(eq_dual, leq_dual)
      #dual_vars <- get_dual_values(solution$y, extract_dual_value, inverse_data[[object@neq_constr]])
      dual_vars <- eq_dual

    }
  } else {
    primal_vars[[as.character(inverse_data[[object@var_id]])]] <- NA_real_
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

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance on the primal and dual residual.
#' @param reltol The relative tolerance on the duality gap.
#' @param abstol The absolute tolerance on the duality gap.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn CPLEX_CONIC Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "CPLEX_CONIC", function(object, data, warm_start, verbose, feastol, reltol, abstol,
                                                    num_iter,
                                                    solver_opts, solver_cache) {
  if (missing(solver_cache)) solver_cache <- new.env(parent=emptyenv())
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
    vtype[data[BOOL_IDX]$bool_vars_idx[[i]]] <- "B"
  }
  #Setting Integer variable types
  for(i in seq_along(data[INT_IDX]$int_vars_idx)){
    vtype[data[BOOL_IDX]$int_vars_idx[[i]]] <- "I"
  }

  #Setting sense of the A matrix
  sense_vec <- rep("E", nrow(Amat))

  #Add inequalities
  leq_start <- dims@zero
  leq_end <- dims@zero + dims@nonpos

  for(j in leq_start:leq_end){
    sense_vec[j + 1] <- "L"
  }

  #Setting Lower bounds of variables
  lb <- rep(-Inf, n_var + total_soc)

  qc <- list()

  soc_start <- leq_start + dims@nonpos
  current_var <- n_var

  for(i in seq_along(dims@soc)){
    for(j in 1:dims@soc[[i]]){
      sense_vec[soc_start + dims@soc[[i]][j]] <- "E"
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

  ## Throw parameter warnings
  if(!all(c(is.null(feastol), is.null(reltol), is.null(abstol)))) {
    warning("Ignoring inapplicable parameter feastol/reltol/abstol for CPLEX.")
  }

  if(is.null(num_iter)) {
      num_iter <- SOLVER_DEFAULT_PARAM$CPLEX$itlim
  }

  #Setting verbosity off
  control <- list(trace = verbose, itlim = num_iter)

  #Setting rest of the parameters
  control[names(solver_opts)] <- solver_opts

  # Solve problem.
  results_dict <- list()

  tryCatch({
    # Define CPLEX problem and solve
    model <- Rcplex::Rcplex_solve_QCP(cvec=cvar, Amat=Amat, bvec=bvec, QC=QC, lb=lb, ub=Inf,
                                      sense=sense_vec, objsense="min", vtype=vtype, control=control)
  }, error = function(e) {
    results_dict$status <- SOLVER_ERROR
  })

  #Changing dualvar to include SOC
  y <- model$extra$lambda
  soc_start <- leq_start + dims@nonpos
  for(i in seq_along(dims@soc)){
    y <- append(y, 0, soc_start)
    soc_start <- soc_start + dims@soc[[i]] + 1
  }
  results_dict$y <- -y
  results_dict$model <- model
  results_dict$eq_dual <- results_dict$y[1:dims@zero]
  results_dict$ineq_dual <- results_dict$y[-(1:dims@zero)]

  return(results_dict)

})

#' An interface for the CVXOPT solver.
#'
setClass("CVXOPT", prototype = list(MIP_CAPABLE = FALSE,
                                    SUPPORTED_CONSTRAINTS = c(supported_constraints(ConicSolver()), "SOC", "PSDConstraint")),
         contains = "ECOS")

CVXOPT <- function() { new("CVXOPT") }

# Map of CVXOPT status to CVXR status.
#' @param solver,object,x A \linkS4class{CVXOPT} object.
#' @param status A status code returned by the solver.
#' @describeIn CVXOPT Converts status returned by the CVXOPT solver to its respective CVXPY status.
setMethod("status_map", "CVXOPT", function(solver, status) {
  if(status == "optimal")
    OPTIMAL
  else if(status %in% c("infeasible", "primal infeasible",
                        "LP relaxation is primal infeasible"))
    INFEASIBLE
  else if(status %in% c("unbounded", "LP relaxation is dual infeasible",
                        "dual infeasible"))
    UNBOUNDED
  else if(status %in% c("solver_error", "unknown", "undefined"))
    SOLVER_ERROR
  else
    stop("CVXOPT status unrecognized: ", status)
})

#' @describeIn CVXOPT Returns the name of the solver.
setMethod("name", "CVXOPT", function(x) { CVXOPT_NAME })

#' @describeIn CVXOPT Imports the solver.

## CVXOPT is not implemented as there is no R package equivalent to cccopt. We should check out cccp, though
setMethod("import_solver", "CVXOPT", function(solver) { requireNamespace("cccp", quietly = TRUE) })

#' @param problem A \linkS4class{Problem} object.
#' @describeIn CVXOPT Can CVXOPT solve the problem?
setMethod("accepts", signature(object = "CVXOPT", problem = "Problem"), function(object, problem) {
  # Can CVXOPT solver the problem?
  # TODO: Check if the matrix is stuffed.
  import_solver(object)
  if(!is_affine(problem@objective@args[[1]]))
    return(FALSE)
  for(constr in problem@constraints) {
    if(!inherits(constr, supported_constraints(object)))
      return(FALSE)
    for(arg in constr@args) {
      if(!is_affine(arg))
        return(FALSE)
    }
  }
  return(TRUE)
})

#' @describeIn CVXOPT Returns a new problem and data for inverting the new solution.
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

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A list containing data necessary for the inversion.
#' @describeIn CVXOPT Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "CVXOPT", solution = "list", inverse_data = "list"), function(object, solution, inverse_data) {
  status <- solution$status
  primal_vars <- list()
  dual_vars <- list()
  if(status %in% SOLUTION_PRESENT){
    opt_val <- solution$value + inverse_data[[OFFSET]]
    primal_vars[[as.character(inverse_data[[object@var_id]])]] <- solution$primal
    eq_dual <- get_dual_values(solution$eq_dual, extract_dual_value, inverse_data[[object@eq_constr]])
    leq_dual <- get_dual_values(solution$ineq_dual, extract_dual_value, inverse_data[[object@neq_constr]])
    eq_dual <- utils::modifyList(eq_dual, leq_dual)
    dual_vars <- eq_dual
    return(Solution(status, opt_val, primal_vars, dual_vars, list()))
  } else {
    return(failure_solution(status))
  }
})

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance on the primal and dual residual.
#' @param reltol The relative tolerance on the duality gap.
#' @param abstol The absolute tolerance on the duality gap.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn CVXOPT Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "CVXOPT", function(object, data, warm_start, verbose, feastol, reltol, abstol,
                                               num_iter, solver_opts, solver_cache) {
  #Tweak parameters
  if(is.null(feastol)) {
    feastol <- SOLVER_DEFAULT_PARAM$CVXOPT$feastol
  }
  if(is.null(reltol)) {
    reltol <- SOLVER_DEFAULT_PARAM$CVXOPT$reltol
  }
  if(is.null(abstol)) {
    abstol <- SOLVER_DEFAULT_PARAM$CVXOPT$abstol
  }
  if(is.null(num_iter)) {
    num_iter <- SOLVER_DEFAULT_PARAM$CVXOPT$max_iters
  }
  param <- cccp::ctrl(maxiters=as.integer(num_iter), abstol=abstol, reltol=reltol,
                      feastol=feastol, trace=as.logical(verbose))
  param$params[names(solver_opts)] <- solver_opts

  G <- as.matrix(data[[G_KEY]])
  h <- as.matrix(data[[H_KEY]])
  nvar <- dim(G)[2]
  dims <- data[[DIMS]]
  zero_dims <- dims@zero
  nonpos_dims <- dims@nonpos
  soc_dims <- dims@soc
  psd_dims <- dims@psd
  clistLength <- (nonpos_dims > 0) + length(soc_dims) + length(psd_dims)

  # For all the constraints except the zero constraint
  clist <- vector(mode="list", length = clistLength)
  clistCounter <- 0
  ghCounter <- 0

  # Deal with non positive constraints
    if(nonpos_dims > 0){
        clistCounter <- clistCounter + 1
        indices  <- seq.int(from = ghCounter + 1, length.out = nonpos_dims)
        clist[[clistCounter]] <- cccp::nnoc(G = G[indices, , drop = FALSE],
                                            h = h[indices, , drop = FALSE])
        ghCounter <- ghCounter + nonpos_dims
    }

  # Deal with SOC constraints
    for(i in soc_dims){
        clistCounter <- clistCounter + 1
        indices  <- seq.int(from = ghCounter + 2, length.out = i - 1)
        clist[[clistCounter]] <- cccp::socc(F = -G[indices, , drop = FALSE],
                                            g = h[indices, , drop = FALSE],
                                            d = -G[ghCounter + 1, , drop = FALSE],
                                            f = h[ghCounter + 1, , drop = FALSE])
        ghCounter <- ghCounter + i
    }

  # Deal with PSD constraints
  for(i in psd_dims){
      Flist <- vector(mode="list", length = nvar+1)
      indices  <- seq.int(from = ghCounter + 1, length.out = i^2)
      currG <- G[indices, , drop = FALSE]
      currh <- h[indices, , drop = FALSE]
      Flist[[1]] <- matrix(currh, nrow = i)
      for(j in seq_len(nvar)){
          Flist[[j+1]] <- matrix(currG[, j, drop = FALSE], nrow = i)
      }
      clistCounter <- clistCounter + 1
      clist[[clistCounter]] <- cccp::psdc(Flist = Flist[-1], F0 = Flist[[1]])
      ghCounter <- ghCounter + i
  }

  if(zero_dims > 0){
    results <- cccp::cccp(q=data[[C_KEY]],
                    A=as.matrix(data[[A_KEY]]),
                    b=as.matrix(data[[B_KEY]]),
                    cList=clist,
                    optctrl=param)
  } else {
    results <- cccp::cccp(q=data[[C_KEY]],
                    cList=clist,
                    optctrl=param)
  }
  solution <- list()
  solution$status <- status_map(object, results$status)
  solution$value <- (results$state[1] + data[[OFFSET]])[[1]]
  solution$primal <- cccp::getx(results)
  solution$eq_dual <- cccp::gety(results)
  solution$ineq_dual <- unlist(cccp::getz(results))
  #solution$ineq_dual <- as.matrix(c(temp[[1]], temp[[2]], temp[[3]], temp[[4]][abs(temp[[4]]) > 1e-8])[1:nrow(G)])

  return(solution)
})


#' An interface for the ECOS BB solver.
#'
#' @name ECOS_BB-class
#' @aliases ECOS_BB
#' @rdname ECOS_BB-class
#' @export
setClass("ECOS_BB", slots = list(MI_SUPPORTED_CONSTRAINTS = "character"),
                    prototype = list(MIP_CAPABLE = TRUE, MI_SUPPORTED_CONSTRAINTS = supported_constraints(ECOS())),
         contains = "ECOS")

#' @rdname ECOS_BB-class
#' @export
ECOS_BB <- function() { new("ECOS_BB") }

#' @param object,x A \linkS4class{ECOS_BB} object.
#' @describeIn ECOS_BB Returns the name of the solver.
setMethod("name", "ECOS_BB", function(x) { ECOS_BB_NAME })

#' @param problem A \linkS4class{Problem} object.
#' @describeIn ECOS_BB Returns a new problem and data for inverting the new solution.
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

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance.
#' @param reltol The relative tolerance.
#' @param abstol The absolute tolerance.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn ECOS_BB Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "ECOS_BB", function(object, data, warm_start, verbose,
                                                feastol,
                                                reltol,
                                                abstol,
                                                num_iter,
                                                solver_opts, solver_cache) {
  if (missing(solver_cache)) solver_cache <- new.env(parent=emptyenv())

  cones <- ECOS.dims_to_solver_dict(data[[ConicSolver()@dims]])

  if(is.null(feastol)) {
      feastol <- SOLVER_DEFAULT_PARAM$ECOS_BB$feastol
  }
  if(is.null(reltol)) {
      reltol <- SOLVER_DEFAULT_PARAM$ECOS_BB$reltol
  }
  if(is.null(abstol)) {
      abstol <- SOLVER_DEFAULT_PARAM$ECOS_BB$abstol
  }
  if(is.null(num_iter)) {
      num_iter <- SOLVER_DEFAULT_PARAM$ECOS_BB$maxit
  }
  num_iter  <- as.integer(num_iter)
  ecos_opts <- ECOSolveR::ecos.control(maxit = num_iter, feastol = feastol, reltol = reltol, abstol = abstol, verbose = as.integer(verbose), mi_max_iters = num_iter)
  ecos_opts[names(solver_opts)] <- solver_opts
  solution <- ECOSolveR::ECOS_csolve(c = data[[C_KEY]], G = data[[G_KEY]], h = data[[H_KEY]], dims = cones, A = data[[A_KEY]], b = data[[B_KEY]],
                                     bool_vars = data[[BOOL_IDX]], int_vars = data[[INT_IDX]], control = ecos_opts)
  return(solution)
})

#' An interface for the GLPK solver.
#'
#' @name GLPK-class
#' @aliases GLPK
#' @rdname GLPK-class
#' @export
setClass("GLPK", prototype = list(MIP_CAPABLE = TRUE,
                                  SUPPORTED_CONSTRAINTS = supported_constraints(ConicSolver())),
         contains = "CVXOPT")

#' @rdname GLPK-class
#' @export
GLPK <- function() { new("GLPK") }

#' @param solver,object,x A \linkS4class{GLPK} object.
#' @param status A status code returned by the solver.
#' @describeIn GLPK Converts status returned by the GLPK solver to its respective CVXPY status.
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

#' @describeIn GLPK Returns the name of the solver.
setMethod("name", "GLPK", function(x) { GLPK_NAME })

#' @describeIn GLPK Imports the solver.
setMethod("import_solver", "GLPK", function(solver) { requireNamespace("Rglpk", quietly = TRUE) })

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A list containing data necessary for the inversion.
#' @describeIn GLPK Returns the solution to the original problem given the inverse_data.
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

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance.
#' @param reltol The relative tolerance.
#' @param abstol The absolute tolerance.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn GLPK Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "GLPK", function(object, data, warm_start, verbose, feastol, reltol, abstol, num_iter, solver_opts, solver_cache) {
  if (missing(solver_cache)) solver_cache <- new.env(parent=emptyenv())
  if(verbose)
    solver_opts$verbose <- verbose
  solver_opts$canonicalize_status <- FALSE

  #Throw warnings if non-default values have been put in
  if(!is.null(feastol)){
    warning("A value has been set for feastol, but the GLPK solver does not accept this parameter. Solver will run without taking this parameter into consideration.")
  }
  if(!is.null(reltol)){
    warning("A value has been set for reltol, but the GLPK solver does not accept this parameter. Solver will run without taking this parameter into consideration.")
  }
  if(!is.null(abstol)){
    warning("A value has been set for abstol, but the GLPK solver does not accept this parameter. Solver will run without taking this parameter into consideration.")
  }
  if(!is.null(num_iter)){
    warning("A value has been set for num_iter, but the GLPK solver does not accept this parameter. Solver will run without taking this parameter into consideration.")
  }

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

#' An interface for the GLPK MI solver.
#'
#' @name GLPK_MI-class
#' @aliases GLPK_MI
#' @rdname GLPK_MI-class
#' @export
setClass("GLPK_MI", prototype = list(MIP_CAPABLE = TRUE,
                                     SUPPORTED_CONSTRAINTS = supported_constraints(ConicSolver()),
                                     MI_SUPPORTED_CONSTRAINTS = supported_constraints(ConicSolver())),
         contains = "GLPK")

#' @rdname GLPK_MI-class
#' @export
GLPK_MI <- function() { new("GLPK_MI") }

# Map of GLPK_MI status to CVXR status.
#' @param solver,object,x A \linkS4class{GLPK_MI} object.
#' @param status A status code returned by the solver.
#' @describeIn GLPK_MI Converts status returned by the GLPK_MI solver to its respective CVXPY status.
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

#' @describeIn GLPK_MI Returns the name of the solver.
setMethod("name", "GLPK_MI", function(x) { GLPK_MI_NAME })

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance.
#' @param reltol The relative tolerance.
#' @param abstol The absolute tolerance.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn GLPK_MI Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "GLPK_MI", function(object, data, warm_start, verbose, feastol, reltol, abstol, num_iter, solver_opts, solver_cache) {
  if (missing(solver_cache)) solver_cache <- new.env(parent=emptyenv())
  if(verbose)
    solver_opts$verbose <- verbose
  solver_opts$canonicalize_status <- FALSE
  if(!is.null(feastol)){
    warning("A value has been set for feastol, but the GLPK solver does not accept this parameter. Solver will run without taking this parameter into consideration.")
  }
  if(!is.null(reltol)){
    warning("A value has been set for reltol, but the GLPK solver does not accept this parameter. Solver will run without taking this parameter into consideration.")
  }
  if(!is.null(abstol)){
    warning("A value has been set for abstol, but the GLPK solver does not accept this parameter. Solver will run without taking this parameter into consideration.")
  }
  if(!is.null(num_iter)){
    warning("A value has been set for num_iter, but the GLPK solver does not accept this parameter. Solver will run without taking this parameter into consideration.")
  }

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

#' An interface for the GUROBI conic solver.
#'
#' @name GUROBI_CONIC-class
#' @aliases GUROBI_CONIC
#' @rdname GUROBI_CONIC-class
#' @export
setClass("GUROBI_CONIC", prototype = list(MIP_CAPABLE = TRUE,
                                          SUPPORTED_CONSTRAINTS = c(supported_constraints(ConicSolver()), "SOC"),
                                          MI_SUPPORTED_CONSTRAINTS = c(supported_constraints(ConicSolver()), "SOC")),
         contains = "SCS")

#' @rdname GUROBI_CONIC-class
#' @export
GUROBI_CONIC <- function() { new("GUROBI_CONIC") }

# Is this one that's used? Should we delete?
# Map of Gurobi status to CVXR status.
# # @describeIn GUROBI_CONIC Converts status returned by the GUROBI solver to its respective CVXPY status.
# setMethod("status_map", "GUROBI_CONIC", function(solver, status) {
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
#     stop("GUROBI status unrecognized: ", status)
# })

#' @param object,x A \linkS4class{GUROBI_CONIC} object.
#' @describeIn GUROBI_CONIC Returns the name of the solver.
setMethod("name", "GUROBI_CONIC", function(x) { GUROBI_NAME })

#' @describeIn GUROBI_CONIC Imports the solver.
setMethod("import_solver", "GUROBI_CONIC", function(solver) { requireNamespace("gurobi", quietly = TRUE) })

# Map of GUROBI status to CVXR status.
#' @param status A status code returned by the solver.
#' @describeIn GUROBI_CONIC Converts status returned by the GUROBI solver to its respective CVXPY status.
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

#' @param problem A \linkS4class{Problem} object.
#' @describeIn GUROBI_CONIC Can GUROBI_CONIC solve the problem?
setMethod("accepts", signature(object = "GUROBI_CONIC", problem = "Problem"), function(object, problem) {
  # TODO: Check if the matrix is stuffed.
  if(!is_affine(problem@objective@args[[1]]))
    return(FALSE)
  for(constr in problem@constraints) {
    if(!inherits(constr, supported_constraints(object)))
      return(FALSE)
    for(arg in constr@args) {
      if(!is_affine(arg))
        return(FALSE)
    }
  }
  return(TRUE)
})

#' @describeIn GUROBI_CONIC Returns a new problem and data for inverting the new solution.
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

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A list containing data necessary for the inversion.
#' @describeIn GUROBI_CONIC Returns the solution to the original problem given the inverse_data.
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

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance.
#' @param reltol The relative tolerance.
#' @param abstol The absolute tolerance.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn GUROBI_CONIC Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "GUROBI_CONIC", function(object, data, warm_start, verbose, feastol, reltol, abstol,
                                                     num_iter, solver_opts, solver_cache) {
  if (missing(solver_cache)) solver_cache <- new.env(parent=emptyenv())
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
  if(!all(c(is.null(reltol), is.null(abstol)))) {
    warning("Ignoring inapplicable parameter reltol/abstol for GUROBI.")
  }
  if(is.null(num_iter)) {
      num_iter <- SOLVER_DEFAULT_PARAM$GUROBI$num_iter
  }
  if (is.null(feastol)) {
      feastol <- SOLVER_DEFAULT_PARAM$GUROBI$FeasibilityTol
  }

  params <- list(OutputFlag = as.numeric(verbose),
                 QCPDual = 1, #equivalent to TRUE
                 IterationLimit = num_iter,
                 FeasibilityTol = feastol,
                 OptimalityTol = feastol)

  params[names(solver_opts)] <- solver_opts

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

