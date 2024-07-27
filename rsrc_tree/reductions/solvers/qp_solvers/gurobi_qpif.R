## CVXPY SOURCE: cvxpy/reductions/solvers/qp_solvers/gurobi_qpif.py
#'
#' An interface for the GUROBI_QP solver.
#'
#' @name GUROBI_QP-class
#' @aliases GUROBI_QP
#' @rdname GUROBI_QP-class
#' @export
setClass("GUROBI_QP", prototype = list(MIP_CAPABLE = TRUE), contains = "QpSolver")

#' @rdname GUROBI_QP-class
#' @export
GUROBI_QP <- function() { new("GUROBI_QP") }

#' @param object,x A \linkS4class{GUROBI_QP} object.
#' @param status A status code returned by the solver.
#' @describeIn GUROBI_QP Converts status returned by the GUROBI solver to its respective CVXPY status.
setMethod("status_map", "GUROBI_QP", function(solver, status) {
  GUROBI_STATUS_MAP <- list(
    "2" = OPTIMAL,
    "3" = INFEASIBLE,
    "5" = UNBOUNDED,
    "4" = INFEASIBLE_OR_UNBOUNDED,
    "6" = INFEASIBLE,
    "7" = SOLVER_ERROR,
    "8" = SOLVER_ERROR,
    "9" = USER_LIMIT,  # Maximum time expired
    # TODO could be anything.
    "10" = SOLVER_ERROR,
    "11" = SOLVER_ERROR,
    "12" = SOLVER_ERROR,
    "13" = OPTIMAL_INACCURATE
  )

  status_string <- GUROBI_STATUS_MAP[[as.character(status)]]
  if (is.null(status_string)) {
      stop("OSQP status unrecognized: ", status)
  } else {
    status_string
  }
})

#' @describeIn GUROBI_QP Returns the name of the solver.
setMethod("name", "GUROBI_QP", function(x) { GUROBI_NAME })

#' @describeIn GUROBI_QP Imports the solver.
setMethod("import_solver", "GUROBI_QP", function(solver) { requireNamespace("gurobi", quietly = TRUE) })


#' @param solution The raw solution returned by the solver.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn GUROBI_QP Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "GUROBI_QP", solution = "list", inverse_data = "InverseData"), function(object, solution, inverse_data){
  model <- solution$model
  solution <- solution$solution
  x_grb <- model$x
  n <- length(x_grb)
  constraints_grb <- model$rhs
  m = length(constraints_grb)

  attr <- list(SOLVE_TIME =  solution$runtime,
               NUM_ITERS = solution$baritercount + solution$itercount)

  status <- status_map(object, solution$status)

  primal_vars <- list()
  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution$objval + inverse_data[[OFFSET]]
    x <- solution$x

    primal_vars[[ names(inverse_data@id_map)[1L] ]] <- x

    #Only add duals if not a MIP
    dual_vars <- list()
    if(!inverse_data@is_mip) {
      if(!is.null(solution$pi)){
        y <- -solution$pi
        dual_vars <- get_dual_values(y, extract_dual_value, inverse_data@sorted_constraints)
      }
    }
  } else {
      primal_vars[[ names(inverse_data@id_map)[1L] ]] <- NA_real_

      dual_vars <- list()
      if(!inverse_data@is_mip) {
        dual_var_ids <- unlist(lapply(inverse_data@sorted_constraints, function(constr) { constr@id }))
        dual_vars <- as.list(rep(NA_real_, length(dual_var_ids)))
        names(dual_vars) <- dual_var_ids
      }

      opt_val <- Inf
      if(status == UNBOUNDED)
        opt_val <- -Inf
  }
  return(Solution(status, opt_val, primal_vars, dual_vars, attr))
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
#' @describeIn GUROBI_QP Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "GUROBI_QP", function(object, data, warm_start, verbose, feastol, reltol, abstol,
                                                  num_iter,
                                                  solver_opts, solver_cache) {
  if (missing(solver_cache)) solver_cache  <- new.env(parent=emptyenv())
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
  vtype <- rep('C', n)

  for(i in seq_along(data[[BOOL_IDX]])){
    vtype[data[[BOOL_IDX]][[i]]] <- 'B' #B for binary
  }

  for(i in seq_along(data[[INT_IDX]])){
    vtype[data[[INT_IDX]][[i]]] <- 'I' #I for integer
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

  # Define objective.

  ####CHECK MATH####
  #Conjecture P is Q matrix and q is c, which is obj for gurobi
  model$Q <- P*.5
  model$obj <- q

  # Throw parameter warnings
  if (!(is.null(reltol) && is.null(abstol))) {
    warning("Ignoring inapplicable parameters reltol/abstol for GUROBI.")
  }

  if (is.null(num_iter)) {
    num_iter <- SOLVER_DEFAULT_PARAM$GUROBI$num_iter
  }

  if (is.null(feastol)) {
      feastol <- SOLVER_DEFAULT_PARAM$GUROBI$FeasibilityTol
  }

  ## Set verbosity and other parameters.
  params <- list(
      OutputFlag = as.numeric(verbose),
      ## TODO: User option to not compute duals.
      QCPDual = 1, #equivalent to TRUE
      IterationLimit = num_iter,
      FeasibilityTol = feastol,
      OptimalityTol = feastol
  )
  params[names(solver_opts)] <- solver_opts

  # Solve problem.
  results_dict <- list()
  tryCatch({
    results_dict$solution <- gurobi::gurobi(model, params)   # Solve.
  }, error = function(e) {   # Error in the solution.
    results_dict$status <- 'SOLVER_ERROR'
  })
  results_dict$model <- model

  results_dict
})

