#'
#' CVXR Package Constants
#'
#' Curvature types as strings
#'
CURV_CONSTANT_KEY = "CONSTANT"
CURV_AFFINE_KEY = "AFFINE"
CURV_CONVEX_KEY = "CONVEX"
CURV_CONCAVE_KEY = "CONCAVE"
CURV_UNKNOWN_KEY = "UNKNOWN"
CURVATURE_STRINGS = c(CURV_CONSTANT_KEY, CURV_AFFINE_KEY, CURV_CONVEX_KEY, CURV_CONCAVE_KEY, CURV_UNKNOWN_KEY)
CURVATURE_NEGATION_MAP <- function(curvature) {
  if(curvature == CURV_CONVEX_KEY)
    CURV_CONCAVE_KEY
  else if(curvature == CURV_CONCAVE_KEY)
    CURV_CONVEX_KEY
  else if(curvature %in% CURVATURE_STRINGS)
    curvature
  else
    stop("Curvature type ", curvature, " not recognized")
}
  
#'
#' Sign types as strings
#'
SIGN_POSITIVE_KEY = "POSITIVE"
SIGN_NEGATIVE_KEY = "NEGATIVE"
SIGN_UNKNOWN_KEY = "UNKNOWN"
SIGN_ZERO_KEY = "ZERO"
SIGN_STRINGS = c(SIGN_POSITIVE_KEY, SIGN_NEGATIVE_KEY, SIGN_UNKNOWN_KEY, SIGN_ZERO_KEY)

#'
#' Monotonicity types as strings
#' 
INCREASING = "INCREASING"
DECREASING = "DECREASING"
SIGNED = "SIGNED"
NONMONOTONIC = "NONMONOTONIC"
MONOTONICITY_STRINGS = c(INCREASING, DECREASING, SIGNED, NONMONOTONIC)

#'
#' Solver capabilities
#'
ECOS.LP_CAPABLE = TRUE
ECOS.SOCP_CAPABLE = TRUE
ECOS.SDP_CAPABLE = FALSE
ECOS.EXP_CAPABLE = TRUE
ECOS.MIP_CAPABLE = FALSE

lp_capable <- function(solver) {
  if(class(solver) == "ECOS")
    return(ECOS.LP_CAPABLE)
  else stop("Unrecognized solver ", name(solver))
}

socp_capable <- function(solver) {
  if(class(solver) == "ECOS")
    return(ECOS.SOCP_CAPABLE)
  else stop("Unrecognized solver ", name(solver))
}

sdp_capable <- function(solver) {
  if(class(solver) == "ECOS")
    return(ECOS.SDP_CAPABLE)
  else stop("Unrecognized solver ", name(solver))
}

exp_capable <- function(solver) {
  if(class(solver) == "ECOS")
    return(ECOS.EXP_CAPABLE)
  else stop("Unrecognized solver ", name(solver))
}

mip_capable <- function(solver) {
  if(class(solver) == "ECOS")
    return(ECOS.LP_CAPABLE)
  else stop("Unrecognized solver ", name(solver))
}

#'
#' Solver exit codes
#'
# EXITCODES from ECOS
# ECOS_OPTIMAL  (0)   Problem solved to optimality
# ECOS_PINF     (1)   Found certificate of primal infeasibility
# ECOS_DINF     (2)   Found certificate of dual infeasibility
# ECOS_INACC_OFFSET (10)  Offset exitflag at inaccurate results
# ECOS_MAXIT    (-1)  Maximum number of iterations reached
# ECOS_NUMERICS (-2)  Search direction unreliable
# ECOS_OUTCONE  (-3)  s or z got outside the cone, numerics?
# ECOS_SIGINT   (-4)  solver interrupted by a signal/ctrl-c
# ECOS_FATAL    (-7)  Unknown problem in solver

# Map of ECOS status to CVXPY status.
ECOS.STATUS_MAP <- function(status) {
  if(status == 0) OPTIMAL
  else if(status == 1) INFEASIBLE
  else if(status == 2) UNBOUNDED
  else if(status == 10) OPTIMAL_INACCURATE
  else if(status == 11) INFEASIBLE_INACCURATE
  else if(status == 12) UNBOUNDED_INACCURATE
  else if(status %in% c(-1, -2, -3, -4, -7)) SOLVER_ERROR
  else stop("ECOS status unrecognized: ", status)
}

status_map <- function(solver, status) {
  if(class(solver) == "ECOS")
    ECOS.STATUS_MAP(status)
  else stop("Unrecognized solver ", name(solver))
}

flatten_list <- function(x) {
  y <- list()
  rapply(x, function(x) y <<- c(y,x))
  y
}

#'
#' Power utilities
#'
is_power2 <- function(num) {
  (round(num) == num) && num > 0 && bitwAnd(num, num - 1) == 0 
}

#'
#' Key utilities
#'
ku_size <- function(key, shape) {
  dims <- c()
  for (i in 1:2) {
    selection <- (1:size(shape)[i])[key[i]]
    size <- size(selection)
    dims <- c(dims, size)
  }
  dims
}