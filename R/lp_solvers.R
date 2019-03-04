CBCSolver <- setClass("CBCSolver", contains = "ReductionSolver")

# Solver capabilities.
setMethod("lp_capable", "CBCSolver", function(object) { TRUE })
setMethod("socp_capable", "CBCSolver", function(object) { FALSE })
setMethod("psd_capable", "CBCSolver", function(object) { FALSE })
setMethod("exp_capable", "CBCSolver", function(object) { FALSE })
setMethod("mip_capable", "CBCSolver", function(object) { TRUE })

# Map of GLPK MIP status to CVXR status.
setMethod("status_map_mip", "CBCSolver", function(object) {
  list(solution = OPTIMAL,
       relaxation_infeasible = INFEASIBLE,
       stopped_on_user_event= SOLVER_ERROR)
})

setMethod("status_map_lp", "CBCSolver", function(object) {
  list(optimal = OPTIMAL,
       primal_infeasible = INFEASIBLE,
       stopped_due_to_errors = SOLVER_ERROR,
       stopped_by_event_handler = SOLVER_ERROR)   # stopped by event handler (virtual int ClpEventHandler::event())
})

setMethod("supported_cut_generators", "CBCSolver", function(object) {
  list(GomoryCuts = "CyCglGomory",
       MIRCuts = "CyCglMixedIntegerRounding",
       MIRCuts2 = "CyCglMixedIntegerRounding2",
       TwoMIRCuts = "CyCglTwomir",
       ResidualCapacityCuts = "CyCglResidualCapacity",
       KnapsackCuts = "CyCglKnapsackCover",
       FlowCoverCuts = "CyCglFlowCover",
       CliqueCuts = "CyCglClique",
       LiftProjectCuts = "CyCglLiftAndProject",
       AllDifferentCuts = "CyCglAllDifferent",
       OddHoleCuts = "CyCglOddHole",
       RedSplitCuts = "CyCglRedSplit",
       LandPCuts = "CyCglLandP",
       PreProcessCuts = "CyCglPreProcess",
       ProbingCuts = "CyCglProbing",
       SimpleRoundingCuts = "CyCglSimpleRounding")
})

# The name of the solver.
setMethod("name", "CBCSolver", function(x) { CBC_NAME })

# Imports the solver.
setMethod("import_solver", "CBCSolver", function(object) {
  if(!requireNamespace("rcbc", quietly = TRUE))
    stop("Required R package rcbc not found. Please install from https://github.com/dirkschumacher/rcbc")
})

# Can CBC solve the problem?
setMethod("accepts", signature(object = "CBCSolver", problem = "Problem"), function(object, problem) {
  # TODO: Check if matrix is stuffed.
  if(!is_affine(problem@objective@args[[1]]))
    return(FALSE)
  for(constr in problem@constraints) {
    if(!class(constr) %in% c("Zero", "NonPos"))
      return(FALSE)
    for(arg in constr@args) {
      if(!is_affine(arg))
        return(FALSE)
    }
  }
  return(TRUE)
})

# Returns a new problem and data for inverting the new solution.
setMethod("perform", signature(object = "CBCSolver", problem = "Problem"), function(object, problem) {
  data <- list()
  inv_data <- list(id(variables(problem)[[1]]))
  names(inv_data) <- object@VAR_ID
  inv_data[OFFSET] <- data[OFFSET][[1]]

  # Order and group constraints.
  eq_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "Zero" })]
  inv_data[object@eq_constr] <- eq_constr
  leq_constr <- problem@constraints[sapply(problem@constraints, function(c) { class(c) == "NonPos" })]
  inv_data[object@neq_constr] <- leq_constr
  return(list(data, inv_data))
})

# Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "CBCSolver", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  status <- solution[STATUS]

  if(status %in% SOLUTION_PRESENT) {
      opt_val <- solution[VALUE]
      primal_vars <- list(solution[PRIMAL])
      names(primal_variables) <- inverse_data[object@var_id]
      eq_dual <- ConicSolver.get_dual_values(solution[EQ_DUAL], inverse_data[EQ_CONSTR])
      leq_dual <- ConicSolver.get_dual_values(solution[INEQ_DUAL], inverse_data[object@neq_constr])
      eq_dual <- modifyList(eq_dual, leq_dual)
      dual_vars <- eq_dual
  } else {
    if(status == INFEASIBLE)
      opt_val <- Inf
    else if(status == UNBOUNDED)
      opt_val <- -Inf
    else
      opt_val <- NA_real_
    primal_vars <- NA
    dual_vars <- NA
  }

  return(Solution(status, opt_val, primal_vars, dual_vars, list()))
})

setMethod("solve", "CBCSolver", function(object, problem, warm_start, verbose, solver_opts) {
  solver <- CBC_OLD()
  inv_data <- perform(object, problem)[[2]]
  objective <- canonical_form(problem@objective)[[1]]
  constraints <- lapply(problem@constraints, function(c) { unlist(canonical_form(c)[[2]]) })
  prob_data <- list(ProblemData())
  names(prob_data) <- name(object)
  sol <- solve(solver, objective, constraints, prob_data, warm_start, verbose, solver_opts)
  invert(object, sol, inv_data)
})
