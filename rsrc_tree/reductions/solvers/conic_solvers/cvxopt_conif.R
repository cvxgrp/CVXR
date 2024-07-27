## CVXPY SOURCE: cvxpy/reductions/solvers/conic_solvers/cvxopt_conif.py

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
