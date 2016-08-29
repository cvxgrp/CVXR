#'
#' Canonical class
#'
setClass("Canonical", contains = "VIRTUAL")
setMethod("canonicalize", "Canonical", function(object) { stop("Unimplemented") })
setMethod("canonical_form", "Canonical", function(object) { canonicalize(object) })
setMethod("variables", "Canonical", function(object) { stop("Unimplemented") })
setMethod("parameters", "Canonical", function(object) { stop("Unimplemented") })
setMethod("constants", "Canonical", function(object) { stop("Unimplemented") })
setMethod("get_data", "Canonical", function(object) { list() })

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
GLPK_NAME = "GLPK"
GLPK_MI_NAME = "GLPK_MI"
CBC_NAME = "CBC"
ECOS_NAME = "ECOS"
ECOS_BB_NAME = "ECOS_BB"
SCS_NAME = "SCS"
GUROBI_NAME = "GUROBI"
MOSEK_NAME = "MOSEK"
LS_NAME = "LS"
SOLVERS <- c(ECOS_NAME)   # TODO: Add more when we implement other solvers

# Parallel (meta) solver
PARALLEL = "parallel"

# Robust CVXOPT LDL KKT solver
ROBUST_KKTSOLVER = "robust"

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
SOLVE_TIME = "solve_time"  # in seconds
SETUP_TIME = "setup_time"  # in seconds
NUM_ITERS = "num_iters"    # number of iterations

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

# Keys for curvature and sign
CONSTANT = "CONSTANT"
AFFINE = "AFFINE"
CONVEX = "CONVEX"
CONCAVE = "CONCAVE"
ZERO = "ZERO"
POSITIVE = "POSITIVE"
NEGATIVE = "NEGATIVE"
UNKNOWN = "UNKNOWN"

SIGN_STRINGS = c(ZERO, POSITIVE, NEGATIVE, UNKNOWN)

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
    return(c(lh_shape[1], rh_shape[2]))
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
  length <- prod(size(t))
  SOCAxis(reshape(sum_expr(list(x, y)), c(length, 1)),
          vstack(list(
              reshape(sub_expr(x, y), c(1, length)),
              reshape(mul_expr(two, t, size(t)), c(1, length))
            ), c(2, length)), 0)
}

gm_constrs <- function(t, x_list, p) {
  if(!is_weight(p)) stop("p must be a valid weight vector")
  w <- dyad_completion(p)

  # TODO: I have no idea what the tree decomposition means
  tree <- decompose(w)
  d <- create_var(size(t))
  d[w] <- t
  
  if(length(x_list) < length(w))
    x_list <- c(x_list, t)
  
  if(length(x_list) != length(w))
    stop("Expected length of x_list to be equal to length of w, but got ", length(x_list), " != ", length(w))
  
  for(i in 1:length(w)) {
    p <- w[i]
    v <- x_list[[i]]
    
    if(p > 0) {
      tmp <- rep(0, length(w))
      tmp[i] <- 1
      d[tmp] <- v
    }
  }
  
  constraints <- list()
  nams <- names(tree)   # Need names to allow for numeric "keys" like in Python
  for(i in 1:length(tree)) {
    elem <- nams[i]
    children <- tree[i]
    if(!(1 %in% elem))
      constraints <- c(constraints, list(gm(d[elem], d[children[1]], d[children[2]])))
  }
  constraints
}

get_num <- function(frac) {
  require(MASS)
  if(!is(frac, "fractions")) stop("frac must be of class fractions")
  fchar <- strsplit(attr(frac, "fracs"), "/", fixed = TRUE)[[1]]
  as.numeric(fchar[1])
}

get_denom <- function(frac) {
  require(MASS)
  if(!is(frac, "fractions")) stop("frac must be of class fractions")
  fchar <- strsplit(attr(frac, "fracs"), "/", fixed = TRUE)[[1]]
  if(length(fchar) == 1)
    return(1)
  else
    as.numeric(fchar[2])
}

pow_high <- function(p, max_denom = 1024, cycles = 10) {
  require(MASS)
  if(p <= 1) stop("Must have p > 1")
  p <- fractions(1/p, cycles = cycles, max.denominator = max_denom)
  if(1/p == as.integer(1/p))
    return(list(as.integer(1/p), c(p, 1-p)))
  list(1/p, c(p, 1-p))
}

pow_mid <- function(p, max_denom = 1024, cycles = 10) {
  require(MASS)
  if(p >= 1 || p <= 0) stop("Must have 0 < p < 1")
  p <- fractions(p, cycles = cycles, max.denominator = max_denom)
  list(p, c(p, 1-p))
}

pow_neg <- function(p, max_denom = 1024, cycles = 10) {
  require(MASS)
  if(p >= 0) stop("must have p < 0")
  p <- fractions(p)
  p <- fractions(p/(p-1), cycles = cycles, max.denominator = max_denom)
  list(p/(p-1), c(p, 1-p))
}

is_power2 <- function(num) {
  (round(num) == num) && num > 0 && bitwAnd(num, num - 1) == 0 
}

is_dyad <- function(frac) {
  require(MASS)
  if((round(frac) == frac) && frac >= 0)
    TRUE
  else if(is(frac, "fractions") && frac >= 0 && is_power2(get_denom(frac)))
    TRUE
  else
    FALSE
}

is_dyad_weight <- function(w) {
  is_weight(w) && all(sapply(w, function(f) { is_dyad(f) }))
}

is_weight <- function(w) {
  if(is.matrix(w) || is.vector(w))
    w <- as.list(w)
  valid_elems <- all(sapply(w, function(v) {
    v >= 0 && ((round(v) == v) || is(v, "fractions"))
  }))
  valid_elems && sum(w) == 1
}

next_pow2 <- function(n) {
  if(n <= 0) return(1)
  n2 <- bitwShiftL(1, bit_length(n))  # TODO: Implement bit_length
  if(bitwShiftR(n2, 1) == n)
    return(n)
  else
    return(n2)
}

get_max_denom <- function(tup) {
  max(sapply(tup, function(f) { get_denom(fractions(f)) }))
}

#'
#' Key utilities
#'
validate_key <- function(key, shape) {
  rows <- shape[1]
  cols <- shape[2]
  
  # Change single indices for vectors into double indices
  if(length(key) != 2) {
    if(rows == 1)
      key <- c(c(1), key)
    else if(cols == 1)
      key <- c(key, c(1))
    else
      stop("Invalid index/slice")
  }
  
  # Change numbers into slices and ensure all slices have a start and stop.
  # key <- format_slice(key[1], shape[1]), format_slice(key[2], shape[2])
  key <- mapply(function(slc, dim) { format_slice(slc, dim) }, slc = key, dim = shape)
  names(key) <- c("rows", "cols")
  key
}

to_int <- function(val) { if(is.na(val)) val else as.integer(val) }

wrap_neg_index <- function(index, dim) {
  if(!is.na(index) && index < 0)
    index <- index %% dim
  index
}

index_to_slice <- function(idx) { c(idx, idx+1) }

ku_size <- function(key, shape) {
  dims <- c()
  for (i in 1:2) {
    selection <- (1:size(shape)[i])[key[i]]
    size <- size(selection)
    dims <- c(dims, size)
  }
  dims
}

index_to_slice <- function(idx) {
  c(idx, idx+1)
}
