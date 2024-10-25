## CVXPY SOURCE: cvxpy/reductions/solvers/conic_solvers/ecos_bb.py


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
