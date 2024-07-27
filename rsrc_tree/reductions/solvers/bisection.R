## CVXPY SOURCE: cvxpy/reductions/solvers/bisection.py

################################################
#                 BISECTION
################################################
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
