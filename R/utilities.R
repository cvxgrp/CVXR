#'
#' CVXR Package Constants
#'
# Constants for operators.
PLUS = "+"
MINUS = "-"
MUL = "*"

# Prefix for default named variables.
VAR_PREFIX = "var"
# Prefix for default named parameters.
PARAM_PREFIX = "param"

# Constraint types
EQ_CONSTR = "=="
INEQ_CONSTR = "<="

# Solver Constants
OPTIMAL = "optimal"
OPTIMAL_INACCURATE = "optimal_inaccurate"
INFEASIBLE = "infeasible"
INFEASIBLE_INACCURATE = "infeasible_inaccurate"
UNBOUNDED = "unbounded"
UNBOUNDED_INACCURATE = "unbounded_inaccurate"
SOLVER_ERROR = "solver_error"
# Statuses that indicate a solution was found.
SOLUTION_PRESENT = c(OPTIMAL, OPTIMAL_INACCURATE)
# Statuses that indicate the problem is infeasible or unbounded.
INF_OR_UNB = c(INFEASIBLE, INFEASIBLE_INACCURATE, UNBOUNDED, UNBOUNDED_INACCURATE)

# Solver names.
CVXOPT_NAME = "CVXOPT"
ECOS_NAME = "ECOS"
ECOS_BB_NAME = "ECOS_BB"
SOLVERS <- c(ECOS_NAME)

# Parallel (meta) solver
PARALLEL = "parallel"

# Map of constraint types
EQ_MAP = 1
LEQ_MAP = 2
SOC_MAP = 3
SOC_EW_MAP = 4
SDP_MAP = 5
EXP_MAP = 6
BOOL_MAP = 7
INT_MAP = 8

# Keys in the dictionary of cone dimensions.
EQ_DIM = "f"
LEQ_DIM = "l"
SOC_DIM = "q"
SDP_DIM = "s"
EXP_DIM = "ep"
# Keys for non-convex constraints.
BOOL_IDS = "bool_ids"
BOOL_IDX = "bool_idx"
INT_IDS = "int_ids"
INT_IDX = "int_idx"

# Keys for results_dict.
STATUS = "status"
VALUE = "value"
OBJ_OFFSET = "obj_offset"
PRIMAL = "primal"
EQ_DUAL = "eq_dual"
INEQ_DUAL = "ineq_dual"

# Keys for problem data dict.
C = "c"
OFFSET = "offset"
A = "A"
B = "b"
G = "G"
H = "h"
F = "F"
DIMS = "dims"
BOOL_IDX = "bool_vars_idx"
INT_IDX = "int_vars_idx"

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
#' Utility functions for shape
#' 
sum_shapes <- function(shapes) {
  rows <- max(sapply(shapes, function(shape) { shape[1] }))
  cols <- max(sapply(shapes, function(shape) { shape[2] }))
  
  # Validate shapes
  for(shape in shapes) {
    if(!all(shape == c(1,1)) && !all(shape == c(rows,cols)))
      stop("Incompatible dimensions")
  }
  c(rows, cols)
}

mul_shapes <- function(lh_shape, rh_shape) {
  if(all(lh_shape == c(1,1)))
    return(rh_shape)
  else if(all(rh_shape == c(1,1)))
    return(lh_shape)
  else {
    if(lh_shape[2] != rh_shape[1])
      stop("Incompatible dimensions")
    return(lh_shape[1], rh_shape[2])
  }
}

#' 
#' Utility functions for sign
#' 
sum_signs <- function(exprs) {
  is_pos <- all(sapply(exprs, function(expr) { is_positive(expr) }))
  is_neg <- all(sapply(exprs, function(expr) { is_negative(expr) }))
  c(is_pos, is_neg)
}

mul_sign <- function(lh_expr, rh_expr) {
  # ZERO * ANYTHING == ZERO
  # POSITIVE * POSITIVE == POSITIVE
  # NEGATIVE * POSITIVE == NEGATIVE
  # NEGATIVE * NEGATIVE == POSITIVE
  is_pos <- (is_zero(lh_expr) || is_zero(rh_expr)) || (is_positive(lh_expr) && is_positive(rh_expr)) || (is_negative(lh_expr) && is_negative(rh_expr))
  is_neg <- (is_zero(lh_expr) || is_zero(rh_expr)) || (is_positive(lh_expr) && is_negative(rh_expr)) || (is_negative(lh_expr) && is_positive(rh_expr))
  c(is_pos, is_neg)
}

#'
#' Utility functions for constraints
#' 
format_elemwise <- function(vars_) {
  spacing <- length(vars_)
  prod_size <- c(spacing * vars_[[1]]$size[1], vars_[[1]]$size[2])
  mat_size <- c(spacing * vars_[[1]]$size[1], vars_[[1]]$size[1])
  
  mat <- lapply(1:spacing, function(i) { get_spacing_matrix(mat_size, spacing, i) })
  terms <- lapply(vars_, function(var) { mul_expr(mat, var, prod_size) })
  list(create_geq(sum_expr(terms)))
}

get_spacing_matrix <- function(size, spacing, offset) {
  require(Matrix)
  col_arr <- 1:size[2]
  row_arr <- spacing * col_arr + offset
  val_arr <- rep(1.0, size[2])
  mat <- sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = size)
  create_const(mat, size, sparse = TRUE)
}

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
gm <- function(t, x, y) {
  two <- create_const(2, c(1, 1))
  SOCElemwise(sum_expr(list(x, y)), list(sub_expr(x, y), mul_expr(two, t, t$size)))  
}

gm_constrs <- function(t, x_list, p) {
  if(!is_weight(p)) stop("p must be a valid weight vector")
  w <- dyad_completion(p)
  
  # TODO: I have no idea what the tree decomposition means
  
  if(length(x_list) < length(w))
    x_list <- c(x_list, t)
  
  if(length(x_list) != length(w))
    stop("Expected length of x_list to be equal to length of w, but got ", length(x_list), " != ", length(w))
  
  # TODO: Iterate through elements in tree
  
  constraints <- list()
  constraints
}

pow_high <- function(p, max_denom = 1024, cycles = 3) {
  require(MASS)
  if(p <= 1) stop("Must have p > 1")
  p <- fractions(1/p, cycles = cycles, max.denominator = max_denom)
  list(1/p, c(p, 1-p))
}

pow_mid <- function(p, max_denom = 1024, cycles = 3) {
  require(MASS)
  if(p >= 1 || p <= 0) stop("Must have 0 < p < 1")
  p <- fractions(p, cycles = cycles, max.denominator = max_denom)
  list(p, c(p, 1-p))
}

pow_neg <- function(p, max_denom = 1024, cycles = 3) {
  require(MASS)
  if(p >= 0) stop("must have p < 0")
  p <- fractions(p/(p-1), cycles = cycles, max.denominator = max_denom)
  list(p/(p-1), c(p, 1-p))
}

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

wrap_neg_index <- function(index, dim) {
  if(!is.na(index) && index < 0)
    index <- index %% dim
  index
}

index_to_slice <- function(idx) {
  c(idx, idx+1)
}
