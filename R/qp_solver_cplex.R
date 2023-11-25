#'
#' An interface for the CPLEX solver.
#'
#' @name CPLEX_QP-class
#' @aliases CPLEX_QP
#' @rdname CPLEX_QP-class
#' @export
setClass("CPLEX_QP", prototype = list(MIP_CAPABLE = TRUE), contains = "QpSolver")

#' @rdname CPLEX_QP-class
#' @export
CPLEX_QP <- function() { new("CPLEX_QP") }

constraint_cplex_infty <- function(v) {
  # Limit values of vector v between +/- infinity as defined in the CPLEX library to be 1e20.
  # https://www.ibm.com/docs/en/icos/22.1.1?topic=keywords-infinity
  CPLEX_INF <- 1e20
  # v[v >= CPLEX_INF] <- CPLEX_INF
  # v[v <= -CPLEX_INF] <- -CPLEX_INF
  return(pmin(pmax(v, -CPLEX_INF), CPLEX_INF))
}

#' @param x,object A \linkS4class{CPLEX_QP} object.
#' @param status A status code returned by the solver.
#' @param default A status string to return if no status code match is found. If \code{default = NA}, this method will return an error when there is no match.
#' @describeIn CPLEX_QP Converts status returned by the CPLEX solver to its respective CVXPY status.
setMethod("status_map", "CPLEX_QP", function(solver, status) {
  CPLEX_STATUS_MAP <- list(
    "1" = OPTIMAL,
    "101" = OPTIMAL,
    "102" = OPTIMAL,
    "3" = INFEASIBLE,
    "4" = INFEASIBLE,
    "22" = INFEASIBLE,
    "103" = INFEASIBLE,
    "2" = UNBOUNDED,
    "21" = UNBOUNDED,
    "118" = UNBOUNDED,
    "4" = UNBOUNDED_INACCURATE,
    "10" = USER_LIMIT
  )

  status_string <- CPLEX_STATUS_MAP[[as.character(status)]]
  if (is.null(status_string)) {
      stop("OSQP status unrecognized: ", status)
  } else {
    status_string
  }
})

#' @describeIn CPLEX_QP Returns the name of the solver.
setMethod("name", "CPLEX_QP", function(x) { CPLEX_NAME })

#' @describeIn CPLEX_QP Imports the solver.
setMethod("import_solver", "CPLEX_QP", function(solver) { requireNamespace("Rcplex", quietly = TRUE) })

#' @param results The raw results returned by the solver.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn CPLEX_QP Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "CPLEX_QP", solution = "list", inverse_data = "InverseData"), function(object, solution, inverse_data){
  model <- solution$model
  attr <- list()

  # TODO: Can't seem to find a way to increase verbosity of cplex. Can't get cputime
  # if("cputime" %in% names(solution))
  #   attr[[SOLVE_TIME]] <- solution$cputime
  # if(inverse_data[[object@IS_MIP]])
  #   attr[[NUM_ITERS]] <- 0
  # else
  #   attr[[NUM_ITERS]] <- as.integer(get_num_barrier_iterations(model$solution$progress))

  status <- status_map(object, solution$model$status, SOLVER_ERROR)

  if(status %in% SOLUTION_PRESENT) {
    # Get objective value.
    opt_val <- model$obj + inverse_data[[OFFSET]]

    # Get solution.
    primal_vars <- list()
    primal_vars[[object@VAR_ID]] <- model$xopt

    # Only add duals if not a MIP.
    dual_vars <- list()
    if(!inverse_data[[object@IS_MIP]]) {
      y <- -as.matrix(model$extra$lambda)   # There's a negative here, should we keep this?
      dual_vars[[object@DUAL_VAR_ID]] <- y
    }

    sol <- Solution(status, opt_val, primal_vars, dual_vars, attr)
  } else
    sol <- failure_solution(status, attr)
  return(sol)
})

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn CPLEX_QP Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "CPLEX_QP", function(object, data, warm_start, verbose, solver_opts, solver_cache = NULL) {

  if(is.null(solver_cache))
    solver_cache  <- new.env(parent = emptyenv())

  P <- Matrix(data[[P_KEY]], byrow = TRUE, sparse = TRUE)
  q <- data[[Q_KEY]]
  A <- Matrix(data[[A_KEY]], byrow = TRUE, sparse = TRUE)   # Equality constraints
  b <- data[[B_KEY]] #Equality constraint RHS
  Fmat <- Matrix(data[[F_KEY]], byrow = TRUE, sparse = TRUE)   # Inequality constraints
  g <- data[[G_KEY]]   # Inequality constraint RHS
  n_var <- data$n_var
  n_eq <- data$n_eq
  n_ineq <- data$n_ineq

  # Constrain values between bounds.
  b <- constrain_cplex_infty(b)
  g <- constrain_cplex_infty(g)

  #In case the b and g variables are empty
  if(prod(dim(b)) == 0 && prod(dim(g)) == 0)
    bvec <- rep(0, n_var)
  else
    bvec <- c(b, g)

  # Create one big constraint matrix with both inequalities and equalities
  Amat <- rbind(A, Fmat)
  if(n_eq + n_ineq == 0){
    #If both number of equalities and inequalities are 0, then constraints dont matter so set equal
    sense_vec <- c(rep("E", n_var))
  } else{
    sense_vec = c(rep("E", n_eq), rep("L", n_ineq))
  }

  #Initializing variable types
  vtype <- rep("C", n_var)

  #Setting Boolean variable types
  for(i in seq_along(data[BOOL_IDX]$bool_vars_idx)){
    vtype[data[BOOL_IDX]$bool_vars_idx[[i]]] <- "B"
  }
  #Setting Integer variable types
  for(i in seq_along(data[INT_IDX]$int_vars_idx)){
    vtype[data[INT_IDX]$int_vars_idx[[i]]] <- "I"
  }

  # Throw parameter warnings
  if(!is.null(feastol) || !is.null(reltol) || !is.null(abstol)) {
      warning("Ignoring inapplicable parameters feastol/reltol/abstol for CPLEX.")
  }

  if (is.null(num_iter)) {
      num_iter <- SOLVER_DEFAULT_PARAM$CPLEX$itlim
  }
  #Setting verbosity off
  control <- list(trace = verbose, itlim = num_iter)

  #Setting rest of the parameters
  control[names(solver_opts)] <- solver_opts

  # Solve problem.
  results_dict <- list()

  #In case A matrix is empty
  if(prod(dim(Amat)) == 0) {
    Amat <- matrix(0, nrow = length(q), ncol = length(q))
  }

  tryCatch({
    # Define CPLEX problem and solve
    model <- Rcplex::Rcplex(cvec=q, Amat=Amat, bvec=bvec, Qmat=P, lb=-Inf, ub=Inf,
                            control=control, objsense="min", sense=sense_vec, vtype=vtype)
    # control parameter would be used to set specific solver arguments. See cran Rcplex documentation
    }, error = function(e) {
      results_dict$status <- SOLVER_ERROR
    }
  )
  results_dict$model <- model
  return(results_dict)
})

