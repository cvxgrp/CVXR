group_constraints <- function(constraints) {
  constr_map <- list()
  for(c in constraints)
    constr_map[[class(c)]] <- c(constr_map[[class(c)]], c)
  return(constr_map)
}

extract_dual_value <- function(result_vec, offset, constraint) {
  value <- result_vec[offset:(offset + size(constraint))]
  if(size(constraint) == 1)
    value <- as.numeric(value)
  offset <- offset + size(constraint)
  return(list(value, offset))
}

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

ReductionSolver <- setClass("ReductionSolver", contains = "Reduction")

# Solver capabilities.
setMethod("mip_capable", function(object) { FALSE })

# Keys for inverse data.
setMethod("var_id", function(object) { "var_id" })
setMethod("eq_constr", function(object) { "eq_constr" })
setMethod("neq_constr", function(object) { "other_constr" })

setMethod("name", "ReductionSolver", function(x) { stop("Unimplemented") })
setMethod("import_solver", "ReductionSolver", function(object) { stop("Unimplemented") })
setMethod("is_installed", "ReductionSolver", function(object) {
  tryCatch(import_solver(object),
           error = function(e) { return(FALSE) },
           finally = return(TRUE)
          )
})

setMethod("solve_via_data", "ReductionSolver", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  stop("Unimplemented")
})

setMethod("solve", "ReductionSolver", function(object, problem, warm_start, verbose, solver_opts) {
  ret <- apply(object, problem)
  data <- ret[[1]]
  inv_data <- ret[[2]]
  solution <- solve_via_data(object, data, warm_start, verbose, solver_opts)
  return(invert(object, inv_data))
})

ConstantSolver <- setClass("ConstantSolver", contains = "ReductionSolver")
setMethod("mip_capable", "ConstantSolver", function(object) { TRUE })
setMethod("accepts", signature(object = "ConstantSolver", problem = "Problem"), function(object, problem) {
  return(length(variables(problem)) == 0)
})

setMethod("apply", signature(object = "ConstantSolver", problem = "Problem"), function(object, problem) {
  return(list(problem, list()))
})

setMethod("invert", signature(object = "ConstantSolver", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  return(solution)
})

setMethod("name", "ConstantSolver", function(x) { return("CONSTANT_SOLVER") })
setMethod("import_solver", "ConstantSolver", function(object) { })
setMethod("is_installed", "ConstantSolver", function(object) { TRUE })
setMethod("solve_via_data", "ConstantSolver", function(object, data, warm_start, verbose, solver_opts, solver_cache = NA) {
  return(solve(object, data, warm_start, verbose, solver_opts))
})

setMethod("solve", "ConstantSolver", function(object, problem, warm_start, verbose, solver_opts) {
  if(all(sapply(problem@constraints, function(c) { !is.na(value(c)) })))
    return(Solution(OPTIMAL, value(problem@objective), list(), list(), list()))
  else
    return(Solution(INFEASIBLE, NA, list(), list(), list()))
})

construct_solving_chain <- function(problem, solver = NA) {
  # Build a reduction chain from a problem to an installed solver.
  if(!is.na(solver)) {
    if(!(solver %in% INSTALLED_SOLVERS))
      stop(paste("The solver", solver, "is not installed"))
    candidates <- list(solver)
  } else
    candidates <- INSTALLED_SOLVERS
  
  reductions <- list()
  # Evaluate parameters and short-circuit the solver if the problem is constant.
  if(length(parameters(problem)) > 0)
    reductions <- c(reductions, EvalParams())
  if(length(variables(problem)) == 0) {
    reductions <- c(reductions, ConstantSolver())
    return(SolvingChain(reductions = reductions))
  }
  if(accepts(Complex2Real(), problem))
    reductions <- c(reductions, Complex2Real())
  
  # Presently, we have but two reduction chains:
  #   1) Qp2SymbolicQp -> QpMatrixStuffing -> list(a QpSolver)
  #   2) Dcp2Cone -> ConeMatrixStuffing -> list(a ConicSolver)
  # Both of these chains require that the problem be DCP.
  if(!is_dcp(problem))
    stop("Problem does not follow DCP rules")
  
  # Both reduction chains exclusively accept minimization problems.
  if(class(problem@objective) == "Maximize")
    reductions <- c(reductions, FlipObjective())
  
  # Attempt to canonicalize the problem to a linearly constrained QP.
  candidate_qp_solvers <- QP_SOLVERS[sapply(QP_SOLVERS, function(s) { s %in% candidates })]
  # Consider only MIQP solvers if problem is integer.
  if(is_mixed_integer(problem))
    candidate_qp_solvers <- candidate_qp_solvers[sapply(candidate_qp_solvers, function(s) { mip_capable(SOLVER_MAP_QP[s]) })]
  if(length(candidate_qp_solvers) > 0 && accepts(Qp2SymbolicQp(), problem)) {
    idx <- sapply(candidate_qp_solvers, function(s) { min(which(QP_SOLVERS == s)) })
    sorted_candidates[order(idx)]
    solver <- sorted_candidates[1]
    solver_instance <- SOLVER_MAP_QP[solver]
    reductions <- c(reductions, list(CvxAttr2Constr(), Qp2SymbolicQp(), QpMatrixStuffing(), solver_instance))
    return(SolvingChain(reductions = reductions))
  }
  
  candidate_conic_solvers <- candidates[sapply(candidates, function(s) { s %in% CONIC_SOLVERS })]
  if(is_mixed_integer(problem)) {
    candidate_conic_solvers <- candidate_conic_solvers[sapply(candidate_conic_solvers, function(s) { mip_capable(SOLVER_MAP_CONIC[s]) })]
    if(length(candidate_conic_solvers) == 0 && length(candidate_qp_solvers) == 0)
      stop(paste("Problem is mixed-integer, but candidate QP/Conic solvers (",
           paste(c(candidate_qp_solvers, candidate_conic_solvers), sep = "", collapse = ",")
           ,") are not MIP-capable.", sep = ""))
  }
  
  if(length(candidate_conic_solvers) == 0)
    stop(paste("Problem could not be reduced to a QP, and no conic solvers exist among candidate solvers (", 
         paste(candidates, collapse = ","), ")", sep = ""))
  
  # Attempt to canonicalize the problem to a cone program.
  # Our choice of solver depends upon w hich atoms are present in the problem.
  # The types of atoms to check for are SOC atoms, PSD atoms, and exponential atoms.
  atoms <- atoms(problem)
  cones <- list()
  if(any(atoms %in% SOC_ATOMS) || any(sapply(problem@constraints, function(c) { class(c) == "SOC" })))
    cones <- c(cones, SOC)
  if(any(atoms %in% EXP_ATOMS) || any(sapply(problem@constraints, function(c) { class(c) == "ExpCone"})))
    cones <- c(cones, ExpCone)
  if(any(atoms %in% PSD_ATOMS) || any(sapply(problem@constraints, function(c) { class(c) == "PSD" }))
                               || any(sapply(variables(problem), function(v) { is_psd(v) || is_nsd(v) })))
    cones <- c(cones, PSD)
  
  # Here, we make use of the observation that canonicalization only
  # increases the number of constraints in our problem.
  has_constr <- length(cones) > 0 || length(problem@constraints) > 0
  
  idx <- sapply(candidate_conic_solvers, function(s) { min(which(CONIC_SOLVERS == s)) })
  sorted_candidates <- candidate_conic_solvers[order(idx)]
  for(solver in sorted_candidates) {
    solver_instance <- SOLVER_MAP_CONIC[solver]
    if(all(cones %in% supported_constraints(solver_instance)) && (has_constr || !requires_constr(solver_instance))) {
      reductions <- c(reductions, list(Dcp2Cone(), CvxAttr2Constr(), ConeMatrixStuffing(), solver_instance))
      return(SolvingChain(reductions = reductions))
    }
  }
  
  stop(paste("Either candidate conic solvers (", 
       paste(candidate_conic_solvers, sep = " ", collapse = ","), ") do not support the cones output by the problem (",
       paste(sapply(cones, function(cone) { name(cone) }), sep = " ", collapse = ","), 
       "), or there are not enough constraints in the problem.", sep = ""))
}

.SolvingChain <- setClass("SolvingChain", contains = "Chain")
SolvingChain <- function(reductions) { .SolvingChain(reductions = reductions) }

setMethod("initialize", "SolvingChain", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)
  last <- .Object@reductions[length(.Object@reductions)]
  if(!is(last, "Solver"))
    stop("Solving chains must terminate with a Solver.")
  .Object@solver <- last
  return(.Object)
})

setMethod("solve", "SolvingChain", function(object, problem, warm_start, verbose, solver_opts) {
  tmp <- apply(object, problem)
  data <- tmp[[1]]
  inverse_data <- tmp[[2]]
  solution <- solve_via_data(object@solver, data, warm_start, verbose, solver_opts)
  return(invert(object, solution, inverse_data))
})

setMethod("solve_via_data", "SolvingChain", function(object, problem, data, warm_start, verbose, solver_opts) {
  return(solve_via_data(object@solver, data, warm_start, verbose, solver_opts, problem@.solver_cache))
})
