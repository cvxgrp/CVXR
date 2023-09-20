# QPSolver requires objectives to be stuffed in the following way.
#'
#' Is the QP objective stuffed?
#'
#' @param objective A \linkS4class{Minimize} or \linkS4class{Maximize} object representing the optimization objective.
#' @return Is the objective a stuffed QP?
is_stuffed_qp_objective <- function(objective) {
  expr <- expr(objective)
  return(inherits(expr, "AddExpression") && length(expr@args) == 2 && inherits(expr@args[[1]], "QuadForm") && inherits(expr@args[[2]], "MulExpression") && is_affine(expr@args[[2]]))
}

#'
#' A QP solver interface.
#'
setClass("QpSolver", representation(IS_MIP = "character"), prototype(IS_MIP = "IS_MIP"), contains = "ReductionSolver")   # Key IS_MIP for internal use only!

#' @param object A \linkS4class{QpSolver} object.
#' @describeIn QpSolver What classes of constraints does the solver support?
setMethod("supported_constraints", "QpSolver", function(object) { c("ZeroConstraint", "NonPosConstraint") })

#' @describeIn QPSolver Can the solver solve problems that do not have constraints?
setMethod("requires_constr", "QpSolver", function(object) { FALSE })

#' @param object A \linkS4class{QpSolver} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn QpSolver Is this a QP problem?
setMethod("accepts", signature(object = "QpSolver", problem = "Problem"), function(object, problem) {
  return(is(problem, "ParamQuadProg") 
          && (mip_capable(object) || !is_mixed_integer(problem))
          && length(convex_attributes(list(problem@x))) == 0
          && (length(problem@constraints) > 0 || !requires_constr(object))
          && all(sapply(problem@constraints, function(c) { class(c) %in% supported_constraints(object) })))
})

QpSolver.prepare_data_and_inv_data <- function(object, problem) {
  data <- list()
  inv_data <- list()
  inv_data[[object@VAR_ID]] <- id(problem@x)
  
  constr_map <- group_constraints(problem@constraints)
  data[[object@DIMS]] <- data[[object@DIMS]]
  inv_data[[object@DIMS]] <- data[[object@DIMS]]
  
  # Add information about integer variables.
  inv_data[[object@IS_MIP]] <- is_mixed_integer(problem)
  
  data[[PARAM_PROB]] <- problem
  return(list(problem, data, inv_data))
}

#' @describeIn QpSolver Constructs a QP problem data stored in a list
setMethod("perform", signature(object = "QpSolver", problem = "Problem"), function(object, problem) {
  # Construct QP problem data stored in a dictionary.
  # The QP has the following form
  #
  #    minimize 1/2 x' P x + q' x
  #    subject to A x = b
  #               F x <= g
  
  tmp_dat <- prepare_data_and_inv_data(object, problem)
  problem <- tmp_dat[[1]]
  data <- tmp_dat[[2]]
  inv_data <- tmp_dat[[3]]
  
  tmp_parm <- apply_parameters(problem)
  P <- tmp_parm[[1]]
  q <- tmp_parm[[2]]
  d <- tmp_parm[[3]]
  AF <- tmp_parm[[4]]
  bg <- tmp_parm[[5]]
  inv_data[[OFFSET]] <- d
  
  # Get number of variables.
  n <- size(problem@x)
  len_eq <- data[[object@DIMS]]@zero
  len_leq <- data[[object@DIMS]]@nonpos
  
  if(len_eq > 0) {
    A <- AF[1:len_eq,]
    b <- -bg[1:len_eq]
  } else {
    A <- Matrix(0, nrow = 0, ncol = n, sparse = TRUE)
    b <- matrix(0, nrow = 0, ncol = 1)
  }
  
  if(len_leq > 0 && len_leq < nrow(AF)) {
    Fmat <- AF[(len_eq + 1):nrow(AF),]
    g <- -bg[(len_eq + 1):nrow(AF)]
  } else {
    Fmat <- Matrix(0, nrow = 0, ncol = n, sparse = TRUE)
    g <- -matrix(0, nrow = 0, ncol = 1)
  }
  
  # Create dictionary with problem data.
  data[[P_KEY]] <- Matrix(P, sparse = TRUE)
  data[[Q_KEY]] <- q
  data[[A_KEY]] <- Matrix(A, sparse = TRUE)
  data[[B_KEY]] <- b
  data[[F_KEY]] <- Matrix(Fmat, sparse = TRUE)
  data[[G_KEY]] <- g
  data[[BOOL_IDX]] <- sapply(problem@x@boolean_idx, function(t) { t[[1]] })
  data[[INT_IDX]] <- sapply(problem@x@integer_idx, function(t) { t[[1]] })
  data$n_var <- n
  data$n_eq <- nrow(A)
  data$n_ineq <- nrow(Fmat)
  
  return(list(data, inv_data))
})

#'
#' An interface for the CPLEX solver.
#'
#' @name CPLEX_QP-class
#' @aliases CPLEX_QP
#' @rdname CPLEX_QP-class
#' @export
setClass("CPLEX_QP", contains = "QpSolver")

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

#' @param x,object,solver A \linkS4class{CPLEX_QP} object.
#' @describeIn CPLEX_QP Can the solver handle mixed-integer programs?
setMethod("mip_capable", "CPLEX_QP", function(solver) { TRUE })

# TODO: Add more!
#' @param status A status code returned by the solver.
#' @param default A status string to return if no status code match is found. If \code{default = NA}, this method will return an error when there is no match.
#' @describeIn CPLEX_QP Converts status returned by the CPLEX solver to its respective CVXPY status.
setMethod("status_map", "CPLEX_QP", function(solver, status, default = NA_character_) {
  if(status %in% c(1, 101, 102))
    OPTIMAL
  else if(status %in% c(3, 22, 4, 103))
    INFEASIBLE
  else if(status %in% c(2, 21, 118))
    UNBOUNDED
  else if(status %in% c(10, 107))
    USER_LIMIT
  else {
    if(is.na(default))
      stop("CPLEX status unrecognized: ", status)
    else
      default
  }
})

#' @describeIn CPLEX_QP Returns the name of the solver.
setMethod("name", "CPLEX_QP", function(x) { CPLEX_NAME })

#' @describeIn CPLEX_QP Imports the solver.
setMethod("import_solver", "CPLEX_QP", function(solver) { requireNamespace("Rcplex", quietly = TRUE) })

#' @param results The raw results returned by the solver.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn CPLEX_QP Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "CPLEX_QP", results = "list", inverse_data = "InverseData"), function(object, results, inverse_data){
  model <- results$model
  attr <- list()

  # TODO: Can't seem to find a way to increase verbosity of cplex. Can't get cputime
  # if("cputime" %in% names(results))
  #   attr[[SOLVE_TIME]] <- results$cputime
  # if(inverse_data[[object@IS_MIP]])
  #   attr[[NUM_ITERS]] <- 0
  # else
  #   attr[[NUM_ITERS]] <- as.integer(get_num_barrier_iterations(model$solution$progress))

  status <- status_map(object, results$model$status)

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
      dual_vars <- list()
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
  require(Rcplex)
  if(is.null(solver_cache))
    solver_cache  <- new.env(parent=emptyenv())
  
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
  if((0 %in% dim(b)) && (0 %in% dim(g)))
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
  if(!all(c(is.null(feastol), is.null(reltol), is.null(abstol)))) {
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
  if(0 %in% dim(Amat)){
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

#'
#' An interface for the GUROBI_QP solver.
#'
#' @name GUROBI_QP-class
#' @aliases GUROBI_QP
#' @rdname GUROBI_QP-class
#' @export
setClass("GUROBI_QP", contains = "QpSolver")

#' @rdname GUROBI_QP-class
#' @export
GUROBI_QP <- function() { new("GUROBI_QP") }

#' @param solver,object,x A \linkS4class{GUROBI_QP} object.
#' @describeIn GUROBI_QP Can the solver handle mixed-integer programs?
setMethod("mip_capable", "GUROBI_QP", function(solver) { TRUE })

#' @param status A status code returned by the solver.
#' @describeIn GUROBI_QP Converts status returned by the GUROBI solver to its respective CVXPY status.
setMethod("status_map", "GUROBI_QP", function(solver, status) {
  if(status == 2 || status == "OPTIMAL")
    OPTIMAL
  else if(status == 3 || status == 6 || status == "INFEASIBLE") # DK: I added the words because the GUROBI solver seems to return the words
    INFEASIBLE
  else if(status == 5 || status == "UNBOUNDED")
    UNBOUNDED
  else if(status == 4 | status == "INF_OR_UNBD")
    INFEASIBLE_OR_UNBOUNDED
  else if(status %in% c(7, 8, 10, 11, 12))
    SOLVER_ERROR   # TODO: Could be anything
  else if(status == 9)   # Maximum time expired.
    USER_LIMIT
  else if(status == 13)
    OPTIMAL_INACCURATE
  else
    stop("GUROBI status unrecognized: ", status)
})

#' @describeIn GUROBI_QP Returns the name of the solver.
setMethod("name", "GUROBI_QP", function(x) { GUROBI_NAME })

#' @describeIn GUROBI_QP Imports the solver.
setMethod("import_solver", "GUROBI_QP", function(solver) { requireNamespace("gurobi", quietly = TRUE) })

##DK: IS THIS FUNCTION NECESSARY ANYMORE WITH invert?
## @DK:  No, not necessary (BN)
## setMethod("alt_invert", "GUROBI_QP", function(object, results, inverse_data) {
##   model <- results$model
##   x_grb <- getVars(model)
##   n <- length(x_grb)
##   constraints_grb <- getConstrs(model)
##   m <- length(constraints_grb)

##   # Start populating attribute dictionary.
##   attr <- list()
##   attr[[SOLVE_TIME]] <- model@Runtime
##   attr[[NUM_ITERS]] <- model@BarIterCount

##   # Map GUROBI statuses back to CVXR statuses.
##   status <- status_map(object, model@Status)

##   if(status %in% SOLUTION_PRESENT) {
##     opt_val <- model@objVal
##     x <- as.matrix(sapply(1:n, function(i) { x_grb[i]@X }))

##     primal_vars <- list()
##     primal_vars[names(inverse_data@id_map)[1]] <- x

##     # Only add duals if not a MIP.
##     dual_vars <- list()
##     if(!inverse_data@is_mip) {
##       y <- -as.matrix(sapply(1:m, function(i) { constraints_grb[i]@Pi }))
##       dual_vars <- get_dual_values(y, extract_dual_value, inverse_data@sorted_constraints)
##     } else {
##       primal_vars <- list()
##       primal_vars[names(inverse_data@id_map)[1]] <- NA_real_

##       dual_vars <- list()
##       if(!inverse_data@is_mip) {
##         dual_var_ids <- sapply(inverse_data@sorted_constraints, function(constr) { constr@id })
##         dual_vars <- as.list(rep(NA_real_, length(dual_var_ids)))
##         names(dual_vars) <- dual_var_ids
##       }

##       opt_val <- Inf
##       if(status == UNBOUNDED)
##         opt_val <- -Inf
##     }
##   }

##   return(Solution(status, opt_val, primal_vars, dual_vars, attr))
## })

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

  # Throw parameter warnings
  if(!all(c(is.null(reltol), is.null(abstol)))) {
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

  attr <- list()
  attr[[SOLVE_TIME]] <- solution$runtime
  attr[[NUM_ITERS]] <- solution$baritercount

  status <- status_map(object, solution$status)

  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution$objval
    x <- solution$x

    primal_vars <- list()
    primal_vars[[names(inverse_data@id_map)[1]]] <- x

    #Only add duals if not a MIP
    dual_vars <- list()
    if(!inverse_data@is_mip) {
      if(!is.null(solution$pi)){
        y <- -solution$pi
        dual_vars <- get_dual_values(y, extract_dual_value, inverse_data@sorted_constraints)
      }
    }
  } else {
      primal_vars <- list()
      primal_vars[names(inverse_data@id_map)[1]] <- NA_real_

      dual_vars <- list()
      if(!inverse_data@is_mip) {
        dual_var_ids <- sapply(inverse_data@sorted_constraints, function(constr) { constr@id })
        dual_vars <- as.list(rep(NA_real_, length(dual_var_ids)))
        names(dual_vars) <- dual_var_ids
      }

      opt_val <- Inf
      if(status == UNBOUNDED)
        opt_val <- -Inf
  }
  return(Solution(status, opt_val, primal_vars, dual_vars, attr))
})

#'
#' QP interface for the OSQP solver.
#'
#' @name OSQP-class
#' @aliases OSQP
#' @rdname OSQP-class
#' @export
setClass("OSQP", contains = "QpSolver")

#' @rdname OSQP-class
#' @export
OSQP <- function() { new("OSQP") }

#' @param solver,object,x A \linkS4class{OSQP} object.
#' @param status A status code returned by the solver.
#' @param default A status string to return if no status code match is found. If \code{default = NA}, this method will return an error when there is no match.
#' @describeIn OSQP Converts status returned by the OSQP solver to its respective CVXR status.
setMethod("status_map", "OSQP", function(solver, status, default = NA_character_) {
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
  else if(status == -6)
    USER_LIMIT
  else if(status == -2 || status == -5 || status == -10)   # -2: Maxiter reached. -5: Interrupted by user. -10: Unsolved.
    SOLVER_ERROR
 else {
    if(is.na(default))
      stop("OSQP status unrecognized: ", status)
   else
     default
 }
})

#' @describeIn OSQP Returns the name of the solver.
setMethod("name", "OSQP", function(x) { OSQP_NAME })

#' @describeIn OSQP Imports the solver.
## setMethod("import_solver", "OSQP", function(solver) { requireNamespace("osqp", quietly = TRUE) })
## Since OSQP is a requirement, this is always TRUE
setMethod("import_solver", "OSQP", function(solver) { TRUE })

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn OSQP Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "OSQP", solution = "list", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  attr <- list()
  attr[[SOLVE_TIME]] <- solution$info$run_time
  attr[[EXTRA_STATS]] <- solution

  # Map OSQP statuses back to CVXR statuses.
  status <- status_map(object, solution$info$status_val, default = SOLVER_ERROR)
  
  ## DWK CHANGE END
  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution$info$obj_val + inverse_data[[OFFSET]]
    
    primal_vars <- list()
    primal_vars[[object@VAR_ID]] <- as.matrix(solution$x)
    
    dual_vars <- list()
    dual_vars[[object@DUAL_VAR_ID]] <- solution$y
    
    attr[[NUM_ITERS]] <- solution$info$iter
    sol <- Solution(status, opt_val, primal_vars, dual_vars, attr)
  } else
    sol <- failure_solution(status)
  return(sol)
})

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn OSQP Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "OSQP", function(object, data, warm_start, verbose, solver_opts, solver_cache = NULL) {
  require(osqp)
  if(is.null(solver_cache))
    solver_cache  <- new.env(parent = emptyenv())
  
  P <- data[[P_KEY]]
  q <- data[[Q_KEY]]
  A <- Matrix(do.call(rbind, list(data[[A_KEY]], data[[F_KEY]])), sparse = TRUE)
  data$Ax <- A
  uA <- c(data[[B_KEY]], data[[G_KEY]])
  data$u <- uA
  lA <- c(data[[B_KEY]], rep(-Inf, length(data[[G_KEY]])))
  data$l <- lA

  # Overwrite defaults of eps_abs = eps_rel = 1e-3, max_iter = 4000.
  if(is.null(solver_opts$eps_abs))
      solver_opts$eps_abs <- 1e-5
  if(is.null(solver_opts$eps_rel))
    solver_opts$eps_rel <- 1e-5
  if(is.null(solver_opts$max_iter))
    solver_opts$max_iter <- 10000

  # Use cached data.
  if(warm_start && !is.null(solver_cache) && name(object) %in% names(solver_cache)) {
    cache <- solver_cache[[name(object)]]
    solver <- cache[[1]]
    old_data <- cache[[2]]
    results <- cache[[3]]
    
    new_args <- list()
    for(key in c("q", "l", "u")) {
      if(any(data[[key]] != old_data[[key]]))
        new_args[[key]] <- data[[key]]
    }
    
    factorizing <- FALSE
    if(any(dim(P) != old_data[[P_KEY]]) || any(P != old_data[[P_KEY]])) {
      P_triu <- Matrix::triu(P)
      new_args$Px <- P_triu
      factorizing <- TRUE
    } 
    if(any(dim(A) != old_data$Ax) || any(A != old_data$Ax)) {
      new_args$Ax <- A
      factorizing <- TRUE
    }
    
    if(length(new_args) > 0)
      do.call(solver$Update, new_args)
    
    # Map OSQP statuses back to CVXR statuses.
    status <- status_map(object, results$info$status_val, default = SOLVER_ERROR)
    if(status == OPTIMAL)
      solver$WarmStart(results$x, results$y)
    
    # Polish if factorizing.
    if(is.null(solver_opts$polish))
      solver_opts$polish <- TRUE
    solver$UpdateSettings(c(list(verbose = verbose), solver_opts))
  } else {
    # Initialize and solve problem.
    if(is.null(solver_opts$polish))
      solver_opts$polish <- TRUE
    solver <- osqp::osqp(P, q, A, lA, uA, c(list(verbose = verbose), solver_opts))
  }
  
  results <- solver$Solve()
  
  if(!is.null(solver_cache)) {
    solver_cache[[name(object)]] <- list(solver, data, results)
  }
  return(list(results, solver_cache))
})
