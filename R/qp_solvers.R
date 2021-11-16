# QPSolver requires objectives to be stuffed in the following way.
#'
#' Is the QP objective stuffed?
#'
#' @param objective A \linkS4class{Minimize} or \linkS4class{Maximize} object representing the optimization objective.
#' @return Is the objective a stuffed QP?
is_stuffed_qp_objective <- function(objective) {
  expr <- expr(objective)
  return(class(expr) == "AddExpression" && length(expr@args) == 2 && class(expr@args[[1]]) == "QuadForm" && class(expr@args[[2]]) == "MulExpression" && is_affine(expr@args[[2]]))
}

#'
#' A QP solver interface.
#'
setClass("QpSolver", contains = "ReductionSolver")

#' @param object A \linkS4class{QpSolver} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn QpSolver Is this a QP problem?
setMethod("accepts", signature(object = "QpSolver", problem = "Problem"), function(object, problem) {
  return(class(problem@objective) == "Minimize" && is_stuffed_qp_objective(problem@objective) && are_args_affine(problem@constraints) &&
           all(sapply(problem@constraints, function(c) { class(c) == "ZeroConstraint" || class(c) == "NonPosConstraint" })))
})

#' @describeIn QpSolver Constructs a QP problem data stored in a list
setMethod("perform", signature(object = "QpSolver", problem = "Problem"), function(object, problem) {
  # Construct QP problem data stored in a dictionary.
  # The QP has the following form
  #    minimize 1/2 x' P x + q' x
  #    subject to A x = b
  #               F x <= g
  inverse_data <- InverseData(problem)

  obj <- problem@objective
  # quadratic part of objective is t(x) %*% P %*% x, but solvers expect 0.5*t(x) %*% P %*% x.
  P <- 2*value(expr(obj)@args[[1]]@args[[2]])
  q <- as.vector(value(expr(obj)@args[[2]]@args[[1]]))

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
      coeff_offset <- ConicSolver.get_coeff_offset(expr(con))
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
      coeff_offset <- ConicSolver.get_coeff_offset(expr(con))
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

  return(list(object, data, inverse_data))
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

#' @param x,object,solver A \linkS4class{CPLEX_QP} object.
#' @describeIn CPLEX_QP Can the solver handle mixed-integer programs?
setMethod("mip_capable", "CPLEX_QP", function(solver) { TRUE })

# TODO: Add more!
#' @param status A status code returned by the solver.
#' @describeIn CPLEX_QP Converts status returned by the CPLEX solver to its respective CVXPY status.
setMethod("status_map", "CPLEX_QP", function(solver, status) {
  if(status %in% c(1, 101, 102))
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

#' @describeIn CPLEX_QP Returns the name of the solver.
setMethod("name", "CPLEX_QP", function(x) { CPLEX_NAME })

#' @describeIn CPLEX_QP Imports the solver.
setMethod("import_solver", "CPLEX_QP", function(solver) { requireNamespace("Rcplex", quietly = TRUE) })

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn CPLEX_QP Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "CPLEX_QP", solution = "list", inverse_data = "InverseData"), function(object, solution, inverse_data){
  model <- solution$model
  attr <- list()

  #Can't seem to find a way to increase verbosity of cplex. Can't get cputime
  #if("cputime" %in% names(solution))
  #  attr[SOLVE_TIME] <- results$cputime
  #attr[NUM_ITERS] <- as.integer(get_num_barrier_iterations(model@solution@progress))

  status <- status_map(object, solution$model$status)

  if(status %in% SOLUTION_PRESENT) {
    # Get objective value.
    opt_val <- model$obj

    # Get solution.
    primal_vars <- list()
    primal_vars[[names(inverse_data@id_map)[1]]] <- model$xopt

    # Only add duals if not a MIP.
    dual_vars <- list()
    if(!inverse_data@is_mip) {
      y <- -as.matrix(model$extra$lambda) #there's a negative here, should we keep this?
      dual_vars <- get_dual_values(y, extract_dual_value, inverse_data@sorted_constraints)
    }
  } else {
    primal_vars <- list()
    dual_vars <- list()
    opt_val <- Inf
    if(status == UNBOUNDED)
      opt_val <- -Inf
  }

  return(Solution(status, opt_val, primal_vars, dual_vars, attr))
})

## Do we need this function?
## @DK: no, we don't (BN)
## setMethod("invert", "CPLEX_QP", function(object, solution, inverse_data) {
##   model <- solution$model
##   attr <- list()
##   if("cputime" %in% names(solution))
##     attr[[SOLVE_TIME]] <- solution$cputime
##   attr[[NUM_ITERS]] <- as.integer(get_num_barrier_iterations(model@solution@progress))

##   status <- status_map(object, get_status(model@solution))

##   if(status %in% SOLUTION_PRESENT) {
##     # Get objective value.
##     opt_val <- get_objective_value(model@solution)

##     # Get solution.
##     x <- as.matrix(get_values(model@solution))
##     primal_vars <- list()
##     primal_vars[names(inverse_data@id_map)[1]] <- x

##     # Only add duals if not a MIP.
##     dual_vars <- list()
##     if(!inverse_data@is_mip) {
##       y <- -as.matrix(get_dual_values(model@solution))
##       dual_vars <- get_dual_values(y, extract_dual_value, inverse_data@sorted_constraints)
##     }
##   } else {
##     primal_vars <- list()
##     primal_vars[names(inverse_data@id_map)[1]] <- NA_real_

##     dual_vars <- list()
##     if(!inverse_data@is_mip) {
##       dual_var_ids <- sapply(inverse_data@sorted_constraints, function(constr) { constr@id })
##       dual_vars <- as.list(rep(NA_real_, length(dual_var_ids)))
##       names(dual_vars) <- dual_var_ids
##     }

##     opt_val <- Inf
##     if(status == UNBOUNDED)
##       opt_val <- -Inf
##   }

##   return(Solution(status, opt_val, primal_vars, dual_vars, attr))
## })

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance on the primal and dual residual.
#' @param reltol The relative tolerance on the duality gap.
#' @param abstol The absolute tolerance on the duality gap.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn CPLEX_QP Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "CPLEX_QP", function(object, data, warm_start, verbose, feastol, reltol, abstol,
                                                 num_iter,
                                                 solver_opts, solver_cache) {
  if (missing(solver_cache)) solver_cache  <- new.env(parent=emptyenv())
  #P <- Matrix(data[[P_KEY]], byrow = TRUE, sparse = TRUE)
  P <- data[[P_KEY]]
  q <- data[[Q_KEY]]
  #A <- Matrix(data[[A_KEY]], byrow = TRUE, sparse = TRUE) #Equality constraints
  A <- data[[A_KEY]]
  b <- data[[B_KEY]] #Equality constraint RHS
  #Fmat <- Matrix(data[[F_KEY]], byrow = TRUE, sparse = TRUE) #Inequality constraints
  Fmat <- data[[F_KEY]]
  g <- data[[G_KEY]] #inequality constraint RHS
  n_var <- data$n_var
  n_eq <- data$n_eq
  n_ineq <- data$n_ineq

  #In case the b and g variables are empty
  if( (0 %in% dim(b)) & (0 %in% dim(g)) ){
    bvec <- rep(0, n_var)
  } else{
    bvec <- c(b, g)
  }

  #Create one big constraint matrix with both inequalities and equalities
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
    #control parameter would be used to set specific solver arguments. See cran Rcplex documentation
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
  else if(status == 3 || status == 6 || status == "INFEASIBLE") #DK: I added the words because the GUROBI solver seems to return the words
    INFEASIBLE
  else if(status == 5 || status == "UNBOUNDED")
    UNBOUNDED
  else if(status == 4 | status == "INF_OR_UNBD")
    INFEASIBLE_INACCURATE
  else if(status %in% c(11,12) || status %in% c("INTERRUPTED", "NUMERIC"))
    SOLVER_ERROR   # TODO: Could be anything
  else if(status %in% c(7,8,9,10,13) || status %in% c("ITERATION_LIMIT", "NODE_LIMIT", "TIME_LIMIT", "SOLUTION_LIMIT", "SUBOPTIMAL"))
    OPTIMAL_INACCURATE   # user limit exceeded (but suboptimal solution may still exist) or known to be suboptimal
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
#' An interface for the OSQP solver.
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
#' @describeIn OSQP Converts status returned by the OSQP solver to its respective CVXPY status.
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

#' @describeIn OSQP Returns the name of the solver.
setMethod("name", "OSQP", function(x) { OSQP_NAME })

#' @describeIn OSQP Imports the solver.
##setMethod("import_solver", "OSQP", function(solver) { requireNamespace("osqp", quietly = TRUE) })
## Since OSQP is a requirement, this is always TRUE
setMethod("import_solver", "OSQP", function(solver) { TRUE })

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn OSQP Returns the solution to the original problem given the inverse_data.
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

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance.
#' @param reltol The relative tolerance.
#' @param abstol The absolute tolerance.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn OSQP Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "OSQP", function(object, data, warm_start, verbose, feastol,
                                             reltol,
                                             abstol,
                                             num_iter,
                                             solver_opts, solver_cache) {
  if (missing(solver_cache)) solver_cache  <- new.env(parent=emptyenv())
  P <- data[[P_KEY]]
  q <- data[[Q_KEY]]
  A <- Matrix(do.call(rbind, list(data[[A_KEY]], data[[F_KEY]])), sparse = TRUE)
  data$full_A <- A
  uA <- c(data[[B_KEY]], data[[G_KEY]])
  data$u <- uA
  lA <- c(data[[B_KEY]], rep(-Inf, length(data[[G_KEY]])))
  data$l <- lA

  if(is.null(feastol)) {
      feastol <- SOLVER_DEFAULT_PARAM$OSQP$eps_prim_inf
  }
  if(is.null(reltol)) {
      reltol <- SOLVER_DEFAULT_PARAM$OSQP$eps_rel
  }
  if(is.null(abstol)) {
      abstol <- SOLVER_DEFAULT_PARAM$OSQP$eps_abs
  }
  if(is.null(num_iter)) {
      num_iter <- SOLVER_DEFAULT_PARAM$OSQP$max_iter
  } else {
      num_iter <- as.integer(num_iter)
  }

  control <- list(max_iter = num_iter,
                  eps_abs = abstol,
                  eps_rel = reltol,
                  eps_prim_inf = feastol,
                  eps_dual_inf = feastol,
                  verbose = verbose)

  control[names(solver_opts)] <- solver_opts

  if(length(solver_cache$cache) > 0 && name(object) %in% names(solver_cache$cache)) {
    # Use cached data.
    cache <- solver_cache$cache[[name(object)]]
    solver <- cache[[1]]
    old_data <- cache[[2]]
    results <- cache[[3]]

    same_pattern <- all(dim(P) == dim(old_data[[P_KEY]])) &&
        all(P@p == old_data[[P_KEY]]@p) &&
        all(P@i == old_data[[P_KEY]]@i) &&
        all(dim(A) == dim(old_data$full_A)) &&
        all(A@p == old_data$full_A@p) &&
        all(A@i == old_data$full_A@i)
  } else
      same_pattern <- FALSE

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
          do.call(solver$Update, new_args)

      ## Map OSQP statuses back to CVXR statuses.
      #status <- status_map(object, results@info@status_val)
      #if(status == OPTIMAL)
      #    warm_start(solver, results@x, results@y)

      ## Polish if factorizing.
      if(is.null(control$polish)){
        control$polish <- factorizing
      }

      #Update Parameters
      solver$UpdateSettings(control)

  } else {
    if(is.null(control$polish))
      control$polish <- TRUE

    # Initialize and solve problem.
    solver <- osqp::osqp(P, q, A, lA, uA, control)
  }

  results <- solver$Solve()
  if (length(solver_cache) == 0L) {
      solver_cache$cache[[name(object)]] <- list(solver, data, results)
  }

  return(results)
})
