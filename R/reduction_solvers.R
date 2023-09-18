######################################################################
#                      Solver utility functions.
######################################################################

#'
#' Stacks the values of the given variables.
#'
#' @param variables a list of variables.
#' @param default value to use when variable value is NA
#' @param byrow logical value indicating whether to fill the matrix by rows. Defaults to FALSE.
stack_vals <- function(variables, default, byrow = FALSE) {
  value <- list()
  for(variable in variables) {
    val <- value(variable)
    if(is.na(val))   # Unknown values.
      mat <- matrix(rep(default, size(variable)), ncol = 1)
    else
      mat <- matrix(val, ncol = 1, byrow = byrow)
    value <- c(value, list(mat))
  }
  # return(do.call("c", value))
  return(do.call("rbind", value))
}

expcone_permutor <- function(n_cones, exp_cone_order) {
  if(n_cones == 0)
    return(c())
  order <- rep(exp_cone_order, n_cones)   # e.g., c(1,0,2, 1,0,2, 1, 0,2, ...)
  offsets <- 3*do.call("c", lapply(seq(n_cones), function(i) { rep(i-1, 3) }))   # e.g., c(0,0,0, 3,3,3, 6,6,6, ...)
  perm <- order + offsets
  return(perm)
}

# #'
# #' Organize the constraints into a dictionary keyed by constraint names.
# #'
# #' @param constraints a list of constraints.
# #' @return A list of constraint types where constr_map[[cone_type]] maps to a list.
# group_constraints <- function(constraints) {
#  constr_map <- list()
#  for(c in constraints)
#    constr_map[[class(c)]] <- c(constr_map[[class(c)]], c)
#  return(constr_map)
# }

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
#' @param parse_func Function handle for the parser.
#' @param constraints A list of the constraints in the problem.
#' @return A map of constraint ID to dual variable value.
get_dual_values <- function(result_vec, parse_func, constraints) {
  dual_vars <- list()
  offset <- 0
  for(constr in constraints) {
    # TODO: Reshape based on dual variable size.
    parsed <- parse_func(result_vec, offset, constr)
    dual_vars[[as.character(id(constr))]] <- parsed[[1]]
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
#' @param feastol The feasible tolerance.
#' @param reltol The relative tolerance.
#' @param abstol The absolute tolerance.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn ReductionSolver Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "ReductionSolver", function(object, data, warm_start, verbose, feastol, reltol, abstol, num_iter, solver_opts, solver_cache) {
  ##if (missing(solver_cache)) solver_cache  <- new.env(parent=emptyenv())
  stop("Unimplemented")
})

#' @param problem A \linkS4class{Problem} object.
#' @describeIn ReductionSolver Solve a problem represented by data returned from apply.
setMethod("reduction_solve", "ReductionSolver", function(object, problem, warm_start, verbose, feastol, reltol, abstol, num_iter, solver_opts) {
  ret <- perform(object, problem)
  object <- ret[[1]]
  data <- ret[[2]]
  inverse_data <- ret[[3]]
  solution <- solve_via_data(object, data, warm_start, verbose, feastol, reltol, abstol, num_iter, solver_opts)
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
#' @param feastol The feasible tolerance.
#' @param reltol The relative tolerance.
#' @param abstol The absolute tolerance.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn ConstantSolver Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "ConstantSolver", function(object, data, warm_start, verbose, feastol, reltol, abstol, num_iter, solver_opts, solver_cache) {
  ## if (missing(solver_cache)) solver_cache  <- new.env(parent=emptyenv())
  return(reduction_solve(object, data, warm_start, verbose, feastol, reltol, abstol, num_iter, solver_opts))
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
#' @return A \linkS4class{SolvingChain} that can be used to solve the problem.
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
  if(any(atoms %in% SOC_ATOMS) || any(sapply(problem@constraints, inherits, what = "SOC")))
    cones <- c(cones, "SOC")
  if(any(atoms %in% EXP_ATOMS) || any(sapply(problem@constraints, inherits, what = "ExpCone")))
    cones <- c(cones, "ExpCone")
  if(any(atoms %in% PSD_ATOMS) || any(sapply(problem@constraints, inherits, what = "PSDConstraint"))
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
#' @param object A \linkS4class{SolvingChain} object.
#' @param chain A \linkS4class{Chain} to prepend.
#' @describeIn SolvingChain Create and return a new SolvingChain by concatenating chain with this instance.
setMethod("prepend", signature(object = "SolvingChain", chain = "Chain"), function(object, chain) {
  SolvingChain(reductions = c(chain@reductions, object@reductions))
})

#' @param problem The problem to solve.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance.
#' @param reltol The relative tolerance.
#' @param abstol The absolute tolerance.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @describeIn SolvingChain Applies each reduction in the chain to the problem, solves it,
#' and then inverts the chain to return a solution of the supplied problem.
setMethod("reduction_solve", signature(object = "SolvingChain", problem = "Problem"), function(object, problem, warm_start, verbose, feastol, reltol, abstol, num_iter, solver_opts) {
  tmp <- perform(object, problem)
  object <- tmp[[1]]
  data <- tmp[[2]]
  inverse_data <- tmp[[3]]
  solution <- solve_via_data(object@solver, data, warm_start, verbose, feastol, reltol, abstol, num_iter, solver_opts)
  return(invert(object, solution, inverse_data))
})

#' @param data Data for the solver.
#' @describeIn SolvingChain Solves the problem using the data output by the an apply invocation.
setMethod("reduction_solve_via_data", "SolvingChain", function(object, problem, data, warm_start, verbose,
                                                               feastol, reltol, abstol, num_iter, solver_opts) {
  return(solve_via_data(object@solver, data, warm_start, verbose,
                        feastol, reltol, abstol, num_iter, solver_opts, problem@.solver_cache))
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

  if(!gp) {
      if (!is_dcp(problem)) {
          err_msg  <- "Problem does not follow DCP rules."
          if (is_dgp(problem)) {
              err_msg  <- paste(err_msg, "However, the problem does follow DGP rules. Consider calling this function with gp = TRUE")
          }
          stop(err_msg)
      }
  } else {
      if(!is_dgp(problem)) {
          err_msg  <- "Problem does not follow DGP rules."
          if (is_dcp(problem)) {
              err_msg  <- paste(err_msg, "However, the problem does follow DCP rules. Consider calling this function with gp = FALSE")
          }
          stop(err_msg)
      }
  }

  # Dcp2Cone and Qp2SymbolicQp require problems to minimize their objectives.
  if(inherits(problem@objective, "Maximize"))
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

################################################
#              BISECTION
###############################################
Bisection.lower_problem <- function(problem) {
  # Evaluates lazy constraints
  Problem(Minimize(0), c(problem@constraints, lapply(problem@.lazy_constraints, function(con) { con() })))
}

Bisection.infeasible <- function(problem, result) {
  # is.null(problem) || problem@status %in% c(INFEASIBLE, INFEASIBLE_INACCURATE)
  is.null(problem) || result$status %in% c(INFEASIBLE, INFEASIBLE_INACCURATE)
}

# TODO: Eliminate the need for reconstructing the problem (via Bisection.lower_problem). 
# Right now, some constraints are lazy, namely the callable constraints, because for 
# certain values of the parameters they become invalidated. Either support lazy 
# constraints that do not need recompilation, or (if possible) find rewritings that 'just work'.
Bisection.find_bisection_interval <- function(problem, t, solver = NULL, low = NA_real_, high = NA_real_, max_iters = 100) {
  # Finds an interval for bisection
  if(is.na(low))
    low <- ifelse(is_nonneg(t), 0, -1)
  if(is.na(high))
    high <- ifelse(is_nonpos(t), 0, 1)
  
  infeasible_low <- is_nonneg(t)
  feasible_high <- is_nonpos(t)
  for(i in seq(max_iters)) {
    if(!feasible_high) {
      value(t) <- high
      lowered <- Bisection.lower_problem(problem)
      result <- solve(lowered, solver = solver)
      if(Bisection.infeasible(lowered, result)) {
        low <- high
        high <- 2*high
        next
      } else if(result$status %in% SOLUTION_PRESENT)
        feasible_high <- TRUE
      else
        stop("Solver failed with status ", result$status)
    }
    
    if(!infeasible_low) {
      value(t) <- low
      lowered <- Bisection.lower_problem(problem)
      result <- solve(lowered, solver = solver)
      if(Bisection.infeasible(lowered, result))
        infeasible_low <- TRUE
      else if(result$status %in% SOLUTION_PRESENT) {
        high <- low
        low <- 2*low
        next
      } else
        stop("Solver failed with status ", result$status)
    }
    
    if(infeasible_low && feasible_high)
      return(c(low, high))
  }
  
  stop("Unable to find suitable interval for bisection; your problem may be unbounded.")
}

Bisection.bisect_int <- function(problem, solver, t, low, high, tighten_lower, tighten_higher, eps = 1e-6, verbose = FALSE, max_iters = 100) {
  # Bisect problem on the parameter t
  verbose_freq <- 5
  soln <- NA_real_
  
  for(i in seq(max_iters)) {
    if(low > high)
      stop("low must be less than or equal to high")
    
    if(!is.na(soln) && (high - low) <= eps) {
      # the previous iteration might have been infeasible, but
      # the tighten* functions might have narrowed the interval
      # to the optimal value in the previous iteration (hence the
      # soln is not None check)
      return(c(soln, low, high))
    }
    
    query_pt <- (low + high) / 2.0
    if(verbose && i %% verbose_freq == 0) {
      print(paste("(iteration ", i, ") lower bound: ", low, sep = ""))
      print(paste("(iteration ", i, ") upper bound: ", high, sep = ""))
      print(paste("(iteration ", i, ") query point: ", query_pt, sep = ""))
    }
    value(t) <- query_pt
    lowered <- Bisection.lowered_problem(problem)
    result <- solve(lowered, solver = solver)
    
    if(Bisection.infeasible(lowered, result)) {
      if(verbose && i %% verbose_freq == 0)
        print(paste("(iteration ", i, ") query was infeasible"))
      low <- tighten_lower(query_pt)
    } else if(result$status %in% SOLUTION_PRESENT) {
      if(verbose && i %% verbose_freq == 0)
        print(paste("(iteration ", i, ") query was feasible with solution ", result$solution, sep = ""))
      soln <- result$solution
      high <- tighten_higher(query_pt)
    } else {
      if(verbose)
        print("Aborting; the solver failed...")
      stop("Solver failed with status ", result$status)
    }
  }
}

#'
#' Bisection on a one-parameter family of DCP problems.
#' 
#' Bisects on a one-parameter family of DCP problems emitted by Dqcp2Dcp.
#' 
#' @param problem A \linkS4class{Problem} emitted by Dqcp2Dcp.
#' @param solver The solver to use for bisection.
#' @param low (Optional) Lower bound for bisection.
#' @param high (Optional) Upper bound for bisection.
#' @param eps Terminate bisection when width of interval is strictly less than eps.
#' @param verbose A logical value indicating whether to print verbose output related to bisection.
#' @param max_iters The maximum number of iterations to run the bisection.
#' @return A \linkS4class{Solution} object.
setMethod("bisect", "Problem", function(problem, solver = NULL, low = NA_real_, high = NA_real_, eps = 1e-6, verbose = FALSE, max_iters = 100, max_iters_interval_search = 100) {
  if(!.hasSlot(problem, ".bisection_data"))
    stop("bisect only accepts problems emitted by Dqcp2Dcp (with .bisection_data slot)")
  tmp <- problem@.bisection_data
  feas_problem <- tmp[[1]]
  t <- tmp[[2]]
  tighten_lower <- tmp[[3]]
  tighten_higher <- tmp[[4]]
  
  if(verbose)
    cat("\n******************************************************",
        "**************************\n",
        "Preparing to bisect problem\n\n", as.character(Bisection.lower_problem(problem)), "\n", sep = "")
  
  lowered_feas <- Bisection.lower_problem(feas_problem)
  result <- solve(lowered_feas, solver = solver)
  if(Bisection.infeasible(lowered_feas, result)) {
    if(verbose)
      print("Problem is infeasible")
    return(failure_solution(INFEASIBLE))
  }
  
  if(is.na(low) || is.na(high)) {
    if(verbose)
      print("Finding interval for bisection...")
    tmp <- Bisection.find_bisection_interval(problem, t, solver, low, high, max_iters_interval_search)
    low <- tmp[[1]]
    high <- tmp[[2]]
  }
  
  if(verbose) {
    print(paste("initial lower bound:", low))
    print(paste("initial upper bound:", high))
  }
  
  tmp <- Bisection.bisect_int(problem, solver, t, low, high, tighten_lower, tighten_higher, eps, verbose, max_iters)
  soln <- tmp[[1]]
  low <- tmp[[2]]
  high <- tmp[[3]]
  soln@opt_val <- (low + high) / 2.0
  if(verbose)
    cat("Bisection completed, with lower bound ", low, " and upper bound ", high,
          "\n******************************************",
          "**************************************\n", sep = "")
  return(soln)
})
