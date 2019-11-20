# Solver utility functions.

#' 
#' Organize the constraints into a dictionary keyed by constraint names.
#' 
#' @param constraints a list of constraints.
#' @return A list of constraint types where constr_map[[cone_type]] maps to a list.
group_constraints <- function(constraints) {
  constr_map <- list()
  for(c in constraints)
    constr_map[[class(c)]] <- c(constr_map[[class(c)]], c)
  return(constr_map)
}

#' 
#' Gets a specified value of a dual variable.
#' 
#' @param result_vec A vector containing the dual variable values.
#' @param offset An offset to get correct index of dual values.
#' @param constraint A list of the constraints in the problem.
#' @return A list of a dual variable value and its offset.
extract_dual_value <- function(result_vec, offset, constraint) {
  value <- result_vec[seq(offset + 1, length.out = size(constraint))]
  if(size(constraint) == 1)
    value <- as.numeric(value)
  offset <- offset + size(constraint)
  return(list(value, offset))
}

#' 
#' Gets the values of the dual variables.
#' 
#' @param result_vec A vector containing the dual variable values.
#' @param offset An offset to get correct index of dual values.
#' @param constraint A list of the constraints in the problem.
#' @return A map of constrain id to dual variable value.
get_dual_values <- function(result_vec, parse_func, constraints) {
  dual_vars <- list()
  offset <- 0
  for(constr in constraints) {
    # TODO: Reshape based on dual variable size.
    parsed <- parse_func(result_vec, offset, constr)
    dual_vars[[as.character(constr@id)]] <- parsed[[1]]
    offset <- parsed[[2]]
  }
  return(dual_vars)
}

#'
#' The ReductionSolver class.
#'
#' @rdname ReductionSolver-class
setClass("ReductionSolver", representation(var_id = "character", eq_constr = "character", neq_constr = "character"),   # Keys for inverse data. Internal use only!
                            prototype(var_id = "var_id", eq_constr = "eq_constr", neq_constr = "other_constr"), contains = "Reduction")

# Solver capabilities.
#' @param solver,object,x A \linkS4class{ReductionSolver} object.
#' @describeIn ReductionSolver Can the solver handle mixed-integer programs?
setMethod("mip_capable", "ReductionSolver", function(solver) { FALSE })

#' @describeIn ReductionSolver Returns the name of the solver
setMethod("name", "ReductionSolver", function(x) { stop("Unimplemented") })

#' @describeIn ReductionSolver Imports the solver
setMethod("import_solver", "ReductionSolver", function(solver) { stop("Unimplemented") })

#' @describeIn ReductionSolver Is the solver installed?
setMethod("is_installed", "ReductionSolver", function(solver) { import_solver(solver) })

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn ReductionSolver Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "ReductionSolver", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
  stop("Unimplemented")
})

#' @param problem A \linkS4class{Problem} object.
#' @describeIn ReductionSolver Solve a problem represented by data returned from apply.
setMethod("reduction_solve", "ReductionSolver", function(object, problem, warm_start, verbose, solver_opts) {
  ret <- perform(object, problem)
  object <- ret[[1]]
  data <- ret[[2]]
  inverse_data <- ret[[3]]
  solution <- solve_via_data(object, data, warm_start, verbose, solver_opts)
  return(invert(object, solution, inverse_data))
})

#'
#' The ConstantSolver class.
#'
ConstantSolver <- setClass("ConstantSolver", contains = "ReductionSolver")

#' @param solver,object,x A \linkS4class{ConstantSolver} object.
#' @describeIn ConstantSolver Can the solver handle mixed-integer programs?
setMethod("mip_capable", "ConstantSolver", function(solver) { TRUE })

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
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn ConstantSolver Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "ConstantSolver", function(object, data, warm_start, verbose, solver_opts, solver_cache = list()) {
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
    return(Solution(INFEASIBLE, NA, list(), list(), list()))
})


#'
#' Build a reduction chain from a problem to an installed solver.
#' 
#' @param problem The problem for which to build a chain.
#' @param candidates A list of candidate solvers.
#' @return A SolvingChain that can be used to solve the problem.
construct_solving_chain <- function(problem, candidates) {
  reductions <- list()
  if(length(parameters(problem)) > 0)
    reductions <- c(reductions, EvalParams())
  if(length(variables(problem)) == 0) {
    reductions <- c(reductions, ConstantSolver())
    return(SolvingChain(reductions = reductions))
  }
  
  # Conclude the chain with one of the following:
  #   1) QpMatrixStuffing -> [a QpSolver]
  #   2) ConeMatrixStuffing -> [a ConicSolver]
  
  # First, attempt to canonicalize the problem to a linearly constrained QP.
  if(length(candidates$qp_solvers) > 0 && accepts(QpMatrixStuffing(), problem)) {
    idx <- sapply(candidates$qp_solvers, function(s) { min(which(QP_SOLVERS == s)) })
    sorted_candidates <- candidates$qp_solvers[order(idx)]
    solver <- sorted_candidates[[1]]
    solver_instance <- SOLVER_MAP_QP[[solver]]
    reductions <- c(reductions, list(QpMatrixStuffing(), solver_instance))
    return(SolvingChain(reductions = reductions))
  }
  
  if(length(candidates$conic_solvers) == 0)
    stop(paste("Problem could not be reduced to a QP, and no conic solvers exist among candidate solvers (", 
               paste(unlist(candidates), collapse = ","), ")", sep = ""))
  
  # Our choice of solver depends upon which atoms are present in the problem.
  # The types of atoms to check for are SOC atoms, PSD atoms, and exponential atoms.
  atoms <- atoms(problem)
  cones <- c()
  if(any(atoms %in% SOC_ATOMS) || any(sapply(problem@constraints, class) == "SOC"))
    cones <- c(cones, "SOC")
  if(any(atoms %in% EXP_ATOMS) || any(sapply(problem@constraints, class) == "ExpCone"))
    cones <- c(cones, "ExpCone")
  if(any(atoms %in% PSD_ATOMS) || any(sapply(problem@constraints, class) == "PSDConstraint")
                               || any(sapply(variables(problem), function(v) { is_psd(v) || is_nsd(v) })))
    cones <- c(cones, "PSDConstraint")
  
  # Here, we make use of the observation that canonicalization only
  # increases the number of constraints in our problem.
  has_constr <- length(cones) > 0 || length(problem@constraints) > 0
  
  idx <- sapply(candidates$conic_solvers, function(s) { min(which(CONIC_SOLVERS == s)) })
  sorted_candidates <- candidates$conic_solvers[order(idx)]
  for(solver in sorted_candidates) {
    solver_instance <- SOLVER_MAP_CONIC[[solver]]
    if(all(cones %in% supported_constraints(solver_instance)) && (has_constr || !requires_constr(solver_instance))) {
      reductions <- c(reductions, list(ConeMatrixStuffing(), solver_instance))
      return(SolvingChain(reductions = reductions))
    }
  }
  
  stop(paste("Either candidate conic solvers (", 
       paste(candidates$conic_solvers, sep = " ", collapse = ","), ") do not support the cones output by the problem (",
       paste(cones, sep = " ", collapse = ","), "), or there are not enough constraints in the problem.", sep = ""))
}

setClassUnion("ReductionSolverORNULL", c("ReductionSolver", "NULL"))

#'
#' The SolvingChain class.
#'
#' This class represents a reduction chain that ends with a solver.
#'
#' @rdname SolvingChain-class
.SolvingChain <- setClass("SolvingChain", representation(solver = "ReductionSolverORNULL"), prototype(solver = NULL), contains = "Chain")
SolvingChain <- function(problem = NULL, reductions = list()) { .SolvingChain(problem = problem, reductions = reductions) }

setMethod("initialize", "SolvingChain", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)
  if(length(.Object@reductions) == 0)
    stop("Solving chains must terminate with a ReductionSolver")
  
  last <- .Object@reductions[[length(.Object@reductions)]]
  if(!is(last, "ReductionSolver"))
    stop("Solving chains must terminate with a ReductionSolver.")
  .Object@solver <- last
  return(.Object)
})

# Create and return a new SolvingChain by concatenating chain with this instance.
#' @describeIn SolvingChain Create and return a new SolvingChain by concatenating chain with this instance.
setMethod("prepend", signature(object = "SolvingChain", chain = "Chain"), function(object, chain) {
  SolvingChain(reductions = c(chain@reductions, object@reductions))
})

#' @param problem The problem to solve.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param solver_opts A list of Solver specific options
#' @describeIn SolvingChain Applies each reduction in the chain to the problem, solves it,
#' and then inverts the chain to return a solution of the supplied problem.
setMethod("reduction_solve", signature(object = "SolvingChain", problem = "Problem"), function(object, problem, warm_start, verbose, solver_opts) {
  tmp <- perform(object, problem)
  object <- tmp[[1]]
  data <- tmp[[2]]
  inverse_data <- tmp[[3]]
  solution <- solve_via_data(object@solver, data, warm_start, verbose, solver_opts)
  return(invert(object, solution, inverse_data))
})

#' @param problem The problem to solve.
#' @param data Data for the solver.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param solver_opts A list of Solver specific options
#' @describeIn SolvingChain Solves the problem using the data output by the an apply invocation.
setMethod("reduction_solve_via_data", "SolvingChain", function(object, problem, data, warm_start, verbose, solver_opts) {
  return(solve_via_data(object@solver, data, warm_start, verbose, solver_opts, problem@.solver_cache))
})

#'
#' Builds a chain that rewrites a problem into an intermediate representation suitable for numeric reductions.
#' 
#' @param problem The problem for which to build a chain.
#' @param candidates A list of candidate solvers.
#' @param gp A logical value indicating whether the problem is a geometric program.
#' @return A \linkS4class{Chain} object that can be used to convert the problem to an intermediate form.
setMethod("construct_intermediate_chain", signature(problem = "Problem", candidates = "list"), function(problem, candidates, gp = FALSE) {
  reductions <- list()
  if(length(variables(problem)) == 0)
    return(Chain(reductions = reductions))
  
  # TODO: Handle boolean constraints.
  if(Complex2Real.accepts(problem))
    reductions <- c(reductions, Complex2Real())
  if(gp)
    reductions <- c(reductions, Dgp2Dcp())
  
  if(!gp && !is_dcp(problem))
    stop("Problem does not follow DCP rules. However, the problem does follow DGP rules. Consider calling this function with gp = TRUE")
  else if(gp && !is_dgp(problem))
    stop("Problem does not follow DGP rules. However, the problem does follow DCP rules. Consider calling this function with gp = FALSE")
  
  # Dcp2Cone and Qp2SymbolicQp require problems to minimize their objectives.
  if(class(problem@objective) == "Maximize")
    reductions <- c(reductions, FlipObjective())
  
  # First, attempt to canonicalize the problem to a linearly constrained QP.
  if(length(candidates$qp_solvers) > 0 && Qp2SymbolicQp.accepts(problem)) {
    reductions <- c(reductions, list(CvxAttr2Constr(), Qp2SymbolicQp()))
    return(Chain(reductions = reductions))
  }
  
  # Canonicalize it to conic problem.
  if(length(candidates$conic_solvers) == 0)
    stop("Problem could not be reduced to a QP, and no conic solvers exist among candidate solvers")
  reductions <- c(reductions, list(Dcp2Cone(), CvxAttr2Constr()))
  return(Chain(reductions = reductions))
})
