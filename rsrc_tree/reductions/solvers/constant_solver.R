## CVXPY SOURCE: cvxpy/reductions/solvers/constant_solver.py
#'
#' The ConstantSolver class.
#'
ConstantSolver <- setClass("ConstantSolver",
                           prototype = list(MIP_CAPABLE = TRUE),
                           contains = "ReductionSolver")

#' @param problem A \linkS4class{Problem} object.
#' @describeIn ConstantSolver Is the solver capable of solving the problem?
setMethod("accepts", signature(object = "ConstantSolver", problem = "Problem"), function(object, problem) {
  return(length(variables(problem)) == 0)
})

#' @describeIn ConstantSolver Returns a list of the ConstantSolver, Problem, and an empty list.
setMethod("perform", signature(object = "ConstantSolver", problem = "Problem"), function(object, problem) {
  return(list(object, problem, list()))
})

#' @param solution A \linkS4class{Solution} object to invert.
#' @param inverse_data A list containing data necessary for the inversion.
#' @describeIn ConstantSolver Returns the solution.
setMethod("invert", signature(object = "ConstantSolver", solution = "Solution", inverse_data = "list"), function(object, solution, inverse_data) {
  return(solution)
})

#' @describeIn ConstantSolver Returns the name of the solver.
setMethod("name", "ConstantSolver", function(x) { return("CONSTANT_SOLVER") })

#' @describeIn ConstantSolver Imports the solver.
setMethod("import_solver", "ConstantSolver", function(solver) { TRUE })

#' @describeIn ConstantSolver Is the solver installed?
setMethod("is_installed", "ConstantSolver", function(solver) { TRUE })


#' @param data Data for the solver.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance.
#' @param reltol The relative tolerance.
#' @param abstol The absolute tolerance.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn ConstantSolver Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "ConstantSolver", function(object, data, warm_start, verbose, solver_opts, solver_cache = NULL) {
  ## if (missing(solver_cache)) solver_cache  <- new.env(parent=emptyenv())
  return(reduction_solve(object, data, warm_start, verbose, solver_opts))
})

#' @param problem A \linkS4class{Problem} object.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param solver_opts A list of Solver specific options
#' @describeIn ConstantSolver Solve the problem and return a \linkS4class{Solution} object.
setMethod("reduction_solve", "ConstantSolver", function(object, problem, warm_start, verbose, solver_opts) {
  if(all(sapply(problem@constraints, function(c) { !is.na(value(c)) })))
    return(Solution(OPTIMAL, value(problem@objective), list(), list(), list()))
  else
    return(Solution(INFEASIBLE, NA_real_, list(), list(), list()))
})
