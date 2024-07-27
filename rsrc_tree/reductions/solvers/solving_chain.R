## CVXPY SOURCE: cvxpy/reductions/solvers/solving_chain.py

##################################################
#                 SOLVING CHAIN
##################################################

SolvingChain.is_lp <- function(object) {
  # Is problem a linear program?
  for(c in object@constraints) {
    if(!(is(c, "EqConstraint") || is(c, "ZeroConstraint")) || is_pwl(c@args[[1]]))
      return(FALSE)
  }

  for(var in variables(object)) {
    if(is_psd(var) || is_nsd(var))
      return(FALSE)
  }

  return(is_dcp(object) && is_pwl(object@objective@args[[1]]))
}

SolvingChain.solve_as_qp <- function(problem, candidates) {
  conic_not_qp_slvs <- c()
  for(s in candidates$conic_solvers) {
    if(!(s %in% candidates$qp_solvers))
      conic_not_qp_slvs <- c(conic_not_qp_slvs, s)
  }

  if(SolvingChain.is_lp(problem) && length(conic_not_qp_slvs) > 0) {
    # OSQP can take many iterations for LPs; use a conic solver instead.
    # GUROBI and CPLEX QP/LP interfaces are more efficient -> Use them instead of conic if applicable.
    return(FALSE)
  }
  return(length(candidates$qp_solvers) > 0 && Qp2SymbolicQp.accepts(problem))
}

#'
#' Builds a chain that rewrites a problem into an intermediate representation suitable for numeric reductions.
#'
#' @param problem The problem for which to build a chain.
#' @param candidates List of candidate solvers divided into qp_solvers and conic_solvers.
#' @param gp If TRUE, the problem is parsed as a Disciplined Geometric Program instead of a Disciplined Convex Program
#' @return A list of Reduction objects that can be used to convert the problem to an intermediate form.
SolvingChain.reductions_for_problem_class <- function(problem, candidates, gp = FALSE, solver_opts = NULL) {
  reductions <- list()

  # TODO: Handle boolean constraints.
  if(Complex2Real.accepts(problem))
    reductions <- c(reductions, list(Complex2Real()))
  if(gp)
    reductions <- c(reductions, list(Dgp2Dcp()))

  if(!gp && !is_dcp(problem)) {
    append <- build_non_disciplined_error_msg(problem, "DCP")
    if(is_dgp(problem))
      append <- paste(append, "However, the problem does follow DGP rules. Consider calling solve() with gp = TRUE", sep = "\n")
    else if(is_dqcp(problem))
      append <- paste(append, "However, the problem does follow DQCP rules. Consider calling solve() with qcp = TRUE", sep = "\n")
    stop(paste("Problem does not follow DCP rules. Specifically:", append, sep = "\n"))
  } else if(gp && !is_dgp(problem)) {
    append <- build_non_disciplined_error_msg(problem, "DGP")
    if(is_dcp(problem))
      append <- paste(append, "However, the problem does follow DCP rules. Consider calling solve() with gp = FALSE", sep = "\n")
    else if(is_dqcp(problem))
      append <- paste(append, "However, the problem does follow DQCP rules. Consider calling solve() with qcp = TRUE", sep = "\n")
    stop("Problem does not follow DGP rules. ", append)
  }

  # Dcp2Cone and Qp2SymbolicQp require problems to minimize their objectives.
  if(inherits(problem@objective, "Maximize"))
    reductions <- c(reductions, list(FlipObjective()))

  if(!is.null(solver_opts) && "use_quad_obj" %in% names(solver_opts))
    use_quad <- solver_opts$use_quad_obj
  else
    use_quad <- TRUE

  if(SolvingChain.solve_as_qp(problem, candidates) && use_quad)
    reductions <- c(reductions, list(CvxAttr2Constr(), Qp2SymbolicQp()))
  else {
    # Canonicalize into a conic problem.
    if(length(candidates$conic_solvers) == 0)
      stop("Problem could not be reduced to a QP, and no conic solvers exist among candidate solvers")
  }

  constr_types <- sapply(problem@constraints, class)
  if("FiniteSet" %in% constr_types)
    reductions <- c(reductions, list(Valinvec2mixedint()))

  return(reductions)
}

#'
#' Build a reduction chain from a problem to an installed solver.
#'
#' Note that if the supplied problem has zero variables, the solver parameter
#' will be ignored.
#'
#' @param problem The problem for which to build a chain.
#' @param candidates A list of candidate solvers.
#' @param gp If TRUE, the problem is parsed as a Disciplined Geometric Program instead of a Disciplined Convex Program.
#' @param enforce_dpp When TRUE, a DPPError will be throw when trying to parse a non-DPP problem (instead of just a warning). Defaults to FALSE.
#' @param ignore_dpp When TRUE, DPP problems will be treated as non-DPP, which may speed up compilation. Defaults to FALSE.
#' @param canon_backend Specifies which backend to use for canonicalization, which can affect compilation time. Defaults to NA, i.e., selecting the default backend ("CPP").
#' @param solver_opts List of additional arguments to pass to the solver.
#' @return A \linkS4class{SolvingChain} that can be used to solve the problem.
construct_solving_chain <- function(problem, candidates, gp = FALSE, enforce_dpp = FALSE, ignore_dpp = FALSE, canon_backend = NA_character_, solver_opts = NULL) {
  if(length(variables(problem)) == 0)
    return(SolvingChain(reductions = list(ConstantSolver())))
  reductions <- SolvingChain.reductions_for_problem_class(problem, candidates, gp, solver_opts)

  # Process DPP status of the problem.
  dpp_context <- ifelse(gp, "dgp", "dcp")
  dpp_error_msg <- paste("You are solving a parameterized problem that is not DPP. ",
                         "Because the problem is not DPP, subsequent solves will not be ",
                         "faster than the first one. For more information, see the ",
                         "documentation on Discplined Parametrized Programming, at",
                         "\thttps://www.cvxpy.org/tutorial/advanced/index.html#disciplined-parametrized-programming", collapse = "\n")

  if(ignore_dpp || !is_dpp(problem, dpp_context)) {
    # No warning for ignore_dpp.
    if(ignore_dpp)
      reductions <- c(list(EvalParams()), reductions)
    else if(!enforce_dpp) {
      warning(dpp_error_msg)
      reductions <- c(list(EvalParams()), reductions)
    } else
      stop(dpp_error_msg)
  } else if(any(sapply(parameters(problem), is_complex)))
    reductions <- c(list(EvalParams()), reductions)
  else {   # Compilation with DPP.
    parms_ <- parameters(problem)
    n_parameters <- ifelse(length(parms_) > 0, sum(parms_, function(param) { prod(dim(param)) }), 0)
    if(n_parameters >= PARAM_THRESHOLD)
      warning("Your problem has too many parameters for efficient DPP compilation. We suggest setting ignore_dpp = TRUE")
  }

  # Conclude with matrix stuffing; choose one of the following paths:
  #   1) QpMatrixStuffing -> [a QpSolver]
  #   2) ConeMatrixStuffing -> [a ConicSolver]
  if(!is.null(solver_opts) && "use_quad_obj" %in% names(solver_opts))
    use_quad <- solver_opts$use_quad_obj
  else
    use_quad <- TRUE

  if(SolvingChain.solve_as_qp(problem, candidates) && use_quad) {
    # Canonicalize as a QP.
    solver <- candidates$qp_solvers[[1]]
    solver_instance <- SOLVER_MAP_QP[[solver]]
    reductions <- c(reductions, list(QpMatrixStuffing(canon_backend = canon_backend), solver_instance))
    return(SolvingChain(reductions = reductions))
  }

  # Canonicalize as a cone program.
  if(length(candidates$conic_solvers) == 0)
    stop(paste("Problem could not be reduced to a QP, and no conic solvers exist among candidate solvers (",
               paste(unlist(candidates), collapse = ","), ")", sep = ""))

  # We use constr_types to infer an incomplete list of cones that the solver will need after canonicalization.
  constr_types <- sapply(problem@constraints, class)
  constr_types <- unique(constr_types)
  ex_cos <- constr_types[constr_types %in% EXOTIC_CONES]   # The way we populate ex_cos will need to change if and when we have atoms that require exotic cones.
  approx_cos <- constr_types[constr_types %in% APPROX_CONES]
  for(co in ex_cos) {
    sim_cos <- EXOTIC_CONES[[co]]   # Get the set of required simple cones.
    if(!(app_cos %in% constr_types))
      constr_types <- c(constr_types, app_cos)
    constr_types <- constr_types[which(constr_types != co)]
  }

  # We now go over individual elementary cones support by CVXR (SOC, ExpCone,
  # NonNegConstraint, ZeroConstraint, PSDConstraint, PowCone3D) and check if
  # they've appeared in constr_types or if the problem has an atom requiring
  # that cone.
  cones <- c()
  atoms <- atoms(problem)
  if("SOC" %in% constr_types || any(atoms %in% SOC_ATOMS))
    cones <- c(cones, "SOC")
  if("ExpCone" %in% constr_types || any(atoms %in% EXP_ATOMS))
    cones <- c(cones, "ExpCone")
  if(any(sapply(c("IneqConstraint", "NonPosConstraint", "NonNegConstraint"), function(t) { t %in% constr_types } )) || any(atoms %in% NONPOS_ATOMS))
    cones <- c(cones, "NonNegConstraint")
  if("EqConstraint" %in% constr_types || "ZeroConstraint" %in% constr_types)
    cones <- c(cones, "ZeroConstraint")
  if("PSDConstraint" %in% constr_types || any(atoms %in% PSD_ATOMS) || any(variables(problem), function(v) { is_psd(v) || is_nsd(v) }))
    cones <- c(cones, "PSDConstraint")
  if("PowCone3D" %in% constr_types) {
    # If we add in atoms that specifically use the 3D power cone (rather than
    # the ND power cone), then we'll need to check for those atoms here as well.
    cones <- c(cones, "PowCone3D")
  }

  # Here, we make use of the observation that canonicalization only
  # increases the number of constraints in our problem.
  has_constr <- length(cones) > 0 || length(problem@constraints) > 0

  for(solver in candidates$conic_solvers) {
    solver_instance <- SOLVER_MAP_CONIC[[solver]]

    # Cones supported for MI problems may differ from non MI.
    if(is_mixed_integer(problem))
      supported_constraints <- mi_supported_constraints(solver_instance)
    else
      supported_constraints <- supported_constraints(solver_instance)

    if(all(cones %in% supported_constraints) && (has_constr || !requires_constr(solver_instance))) {
      if(length(ex_cost) > 0)
        reductions <- c(reductions, list(Exotic2Common()))
      if("RelEntrConeQuad" %in% approx_cos || "OpRelEntrConeQuad" %in% approx_cos)
        reductions <- c(reductions, list(QuadApprox()))

      # Should the objective be canonicalized to a quadratic?
      if(!is.null(solver_opts) && "use_quad_obj" %in% names(solver_opts))
        use_quad_obj <- solver_opts$use_quad_obj
      else
        use_quad_obj <- TRUE
      quad_obj <- use_quad_obj && supports_quad_obj(solver_instance) && has_quadratic_term(problem@objective@expr)
      reductions <- c(reductions, list(Dcp2Cone(quad_obj = quad_obj), CvxAttr2Constr(), ConeMatrixStuffing(quad_obj = quad_obj, canon_backend = canon_backend), solver_instance))
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

## Add Solving Chain to SolvingChainOrNULL
setIs("SolvingChain", "SolvingChainORNULL")

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

#' @param data Data for the solver.
#' @describeIn SolvingChain Solves the problem using the data output by the an apply invocation.
setMethod("reduction_solve_via_data", "SolvingChain", function(object, problem, data, warm_start, verbose, solver_opts = list()) {
  return(solve_via_data(object@solver, data, warm_start, verbose, solver_opts, problem@.solver_cache))
})

