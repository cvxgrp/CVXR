## CVXPY SOURCE: None
## This is R-specific

##################################
#                                #
# Solver utilities and constants #
#                                #
##################################

# Solver Constants
OPTIMAL <- "optimal"
OPTIMAL_INACCURATE <- "optimal_inaccurate"
INFEASIBLE <- "infeasible"
INFEASIBLE_INACCURATE <- "infeasible_inaccurate"
UNBOUNDED <- "unbounded"
UNBOUNDED_INACCURATE <- "unbounded_inaccurate"
INFEASIBLE_OR_UNBOUNDED <- "infeasible_or_unbounded"
USER_LIMIT <- "user_limit"
SOLVER_ERROR <- "solver_error"

# Statuses that indicate a solution was found.
SOLUTION_PRESENT <- c(OPTIMAL, OPTIMAL_INACCURATE, USER_LIMIT)

# Statuses that indicate the problem is infeasible or unbounded.
INF_OR_UNB <- c(INFEASIBLE, INFEASIBLE_INACCURATE, UNBOUNDED, UNBOUNDED_INACCURATE, INFEASIBLE_OR_UNBOUNDED)

# Statuses that indicate an inaccurate solution.
INACCURATE <- c(OPTIMAL_INACCURATE, INFEASIBLE_INACCURATE, UNBOUNDED_INACCURATE, USER_LIMIT)

# Statuses that indicate an error.
ERROR <- c(SOLVER_ERROR)

## Codes from lpSolveAPI solver (partial at the moment)
DEGENERATE <- "degenerate"
NUMERICAL_FAILURE <- "numerical_failure"
TIMEOUT <- "timeout"
BB_FAILED <- "branch_and_bound_failure"

## Codes from GLPK (partial)
UNDEFINED <- "undefined"

# Solver names.
CBC_NAME <- "CBC"
CLARABEL_NAME <- "CLARABEL"
PIQP_NAME <- "PIQP"
COPT_NAME <- "COPT"
CPLEX_NAME <- "CPLEX"
CVXOPT_NAME <- "CVXOPT"
DIFFCP_NAME <- "DIFFCP"
ECOS_NAME <- "ECOS"
ECOS_BB_NAME <- "ECOS_BB"
GLOP_NAME <- "GLOP"
GLPK_NAME <- "GLPK"
GLPK_MI_NAME <- "GLPK_MI"
GUROBI_NAME <- "GUROBI"
MOSEK_NAME <- "MOSEK"
NAG_NAME <- "NAG"
OSQP_NAME <- "OSQP"
PDLP_NAME <- "PDLP"
PROXQP_NAME <- "PROXQP"
SCIP_NAME <- "SCIP"
SCS_NAME <- "SCS"
SDPA_NAME <- "SDPA"
XPRESS_NAME <- "XPRESS"
SOLVER_NAMES <- c(CLARABEL_NAME, ECOS_NAME, CVXOPT_NAME, GLOP_NAME, GLPK_NAME,
                  GLPK_MI_NAME, SCS_NAME, SDPA_NAME, GUROBI_NAME, OSQP_NAME,
                  CPLEX_NAME, MOSEK_NAME, CBC_NAME, COPT_NAME, XPRESS_NAME,
                  PROXQP_NAME, NAG_NAME, PDLP_NAME, SCIP_NAME)

# Xpress-specific items.
XPRESS_IIS = "XPRESS_IIS"
XPRESS_TROW = "XPRESS_TROW"

# Parallel (meta) solver
PARALLEL <- "parallel"

# Robust CVXOPT LDL KKT solverg
ROBUST_KKTSOLVER <- "robust"

# Solver option defaults
SOLVER_DEFAULT_PARAM <- list(
  OSQP = list(max_iter = 10000, eps_abs = 1e-5, eps_rel = 1e-5, eps_prim_inf = 1e-4),
  PIQP = list(),  ## same as piqp itself
  ECOS = list(maxit = 100, abstol = 1e-8, reltol = 1e-8, feastol = 1e-8),
  ECOS_BB = list(maxit = 1000, abstol = 1e-6, reltol = 1e-3, feastol = 1e-6),
  ## Until cccp fixes the bug I reported, we set the tolerances as below
  CVXOPT = list(max_iters = 100, abstol = 1e-6, reltol = 1e-6, feastol = 1e-6, refinement = 1L, kktsolver = "chol"),
  SCS = list(max_iters = 2500, eps_rel = 1e-4, eps_abs = 1e-4, eps_infeas = 1e-7),
  CPLEX = list(itlim = 10000),
  MOSEK = list(num_iter = 10000),
  GUROBI = list(num_iter = 10000, FeasibilityTol = 1e-6)
)

# Keys for results_dict.
STATUS <- "status"
VALUE <- "value"
OBJ_OFFSET <- "obj_offset"
PRIMAL <- "primal"
EQ_DUAL <- "eq_dual"
INEQ_DUAL <- "ineq_dual"
SOLVER_NAME <- "solver"
SOLVE_TIME <- "solve_time"  # in seconds
SETUP_TIME <- "setup_time"  # in seconds
NUM_ITERS <- "num_iters"    # number of iterations
EXTRA_STATS <- "extra_stats"   # extra solver-specific statistics

#' Update solver options using defaults and return list
#' @param solver_opts the solver options named list
#' @param the defaults the named list to reconcile against
#' @return list of solver options
reconcile_solver_options <- function(solver_opts, defaults) {
  default_opts <- names(defaults)
  for (x in default_opts) {
    if (is.null(solver_opts[[x]])) solver_opts[[x]] <- defaults[[x]]
  }
  solver_opts
}

