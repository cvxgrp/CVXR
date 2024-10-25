## CVXPY SOURCE: cvxpy/reductions/solvers/conic_solvers/glpk_conif.py

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
