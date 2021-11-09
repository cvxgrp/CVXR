##########################
#                        #
# CVXR Package Constants #
#                        #
##########################
# Constants for operators.
PLUS = "+"
MINUS = "-"
MUL = "*"

# Prefix for default named variables.
VAR_PREFIX = "var"
# Prefix for default named parameters.
PARAM_PREFIX = "param"

# Constraint types.
EQ_CONSTR = "=="
INEQ_CONSTR = "<="

# Atom groups.
SOC_ATOMS = c("GeoMean", "Pnorm", "QuadForm", "QuadOverLin", "Power")
EXP_ATOMS = c("LogSumExp", "LogDet", "Entr", "Exp", "KLDiv", "Log", "Log1p", "Logistic")
PSD_ATOMS = c("LambdaMax", "LambdaSumLargest", "LogDet", "MatrixFrac", "NormNuc", "SigmaMax")

# Solver Constants
OPTIMAL = "optimal"
OPTIMAL_INACCURATE = "optimal_inaccurate"
INFEASIBLE = "infeasible"
INFEASIBLE_INACCURATE = "infeasible_inaccurate"
UNBOUNDED = "unbounded"
UNBOUNDED_INACCURATE = "unbounded_inaccurate"
USER_LIMIT <- "user_limit"
SOLVER_ERROR = "solver_error"
# Statuses that indicate a solution was found.
SOLUTION_PRESENT = c(OPTIMAL, OPTIMAL_INACCURATE)
# Statuses that indicate the problem is infeasible or unbounded.
INF_OR_UNB = c(INFEASIBLE, INFEASIBLE_INACCURATE, UNBOUNDED, UNBOUNDED_INACCURATE)
# Statuses that indicate an error.
ERROR <- c(USER_LIMIT, SOLVER_ERROR)

## Codes from lpSolveAPI solver (partial at the moment)
DEGENERATE = "degenerate"
NUMERICAL_FAILURE = "numerical_failure"
TIMEOUT = "timeout"
BB_FAILED = "branch_and_bound_failure"

## Codes from GLPK (partial)
UNDEFINED = "undefined"

# Solver names.
CBC_NAME = "CBC"
CPLEX_NAME = "CPLEX"
CVXOPT_NAME = "CVXOPT"
ECOS_NAME = "ECOS"
ECOS_BB_NAME = "ECOS_BB"
GLPK_NAME = "GLPK"
GLPK_MI_NAME = "GLPK_MI"
GUROBI_NAME = "GUROBI"
JULIA_OPT_NAME = "JULIA_OPT"
MOSEK_NAME = "MOSEK"
OSQP_NAME = "OSQP"
SCS_NAME = "SCS"
SUPER_SCS_NAME = "SUPER_SCS"
XPRESS_NAME = "XPRESS"
# SOLVERS_NAME <- c(ECOS_NAME, ECOS_BB_NAME, SCS_NAME, LPSOLVE_NAME, GLPK_NAME, MOSEK_NAME, GUROBI_NAME)   # TODO: Add more when we implement other solvers

# Solver option defaults
SOLVER_DEFAULT_PARAM <- list(
    OSQP = list(max_iter = 10000, eps_abs = 1e-5, eps_rel = 1e-5, eps_prim_inf = 1e-4),
    ECOS = list(maxit = 100, abstol = 1e-8, reltol = 1e-8, feastol = 1e-8),
    ECOS_BB = list(maxit = 1000, abstol = 1e-6, reltol = 1e-3, feastol = 1e-6),
    ## Until cccp fixes the bug I reported, we set the tolerances as below
    CVXOPT = list(max_iters = 100, abstol = 1e-6, reltol = 1e-6, feastol = 1e-6, refinement = 1L, kktsolver = "chol"),
    SCS = list(max_iters = 2500, eps_rel = 1e-4, eps_abs = 1e-4, eps_infeas = 1e-7),
    CPLEX = list(itlim = 10000),
    MOSEK = list(num_iter = 10000),
    GUROBI = list(num_iter = 10000, FeasibilityTol = 1e-6)
)

# Xpress-specific items.
XPRESS_IIS = "XPRESS_IIS"
XPRESS_TROW = "XPRESS_TROW"

# Parallel (meta) solver
PARALLEL = "parallel"

# Robust CVXOPT LDL KKT solverg
ROBUST_KKTSOLVER = "robust"

# Map of constraint types
EQ_MAP = "1"
LEQ_MAP = "2"
SOC_MAP = "3"
SOC_EW_MAP = "4"
PSD_MAP = "5"
EXP_MAP = "6"
BOOL_MAP = "7"
INT_MAP = "8"

# Keys in the dictionary of cone dimensions.
##SCS3.0 changes this to "z"
##EQ_DIM = "f"
EQ_DIM = "z"
LEQ_DIM = "l"
SOC_DIM = "q"
PSD_DIM = "s"
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
SOLVER_NAME = "solver"
SOLVE_TIME = "solve_time"  # in seconds
SETUP_TIME = "setup_time"  # in seconds
NUM_ITERS = "num_iters"    # number of iterations

# Keys for problem data dict.
C_KEY = "c"
OFFSET = "offset"
P_KEY = "P"
Q_KEY = "q"
A_KEY = "A"
B_KEY = "b"
G_KEY = "G"
H_KEY = "h"
F_KEY = "F"
DIMS = "dims"
BOOL_IDX = "bool_vars_idx"
INT_IDX = "int_vars_idx"

# Keys for curvature and sign
CONSTANT = "CONSTANT"
AFFINE = "AFFINE"
CONVEX = "CONVEX"
CONCAVE = "CONCAVE"
ZERO = "ZERO"
NONNEG = "NONNEGATIVE"
NONPOS = "NONPOSITIVE"
UNKNOWN = "UNKNOWN"

# Keys for log-log curvature
LOG_LOG_CONSTANT = "LOG_LOG_CONSTANT"
LOG_LOG_AFFINE = "LOG_LOG_AFFINE"
LOG_LOG_CONVEX = "LOG_LOG_CONVEX"
LOG_LOG_CONCAVE = "LOG_LOG_CONCAVE"

SIGN_STRINGS = c(ZERO, NONNEG, NONPOS, UNKNOWN)

# Numerical tolerances.
EIGVAL_TOL = 1e-10
PSD_NSD_PROJECTION_TOL = 1e-8
GENERAL_PROJECTION_TOL = 1e-10
SPARSE_PROJECTION_TOL = 1e-10

SolveResult <- list(SolveResult = list("opt_value", "status", "primal_values", "dual_values"))

apply_with_keepdims <- function(x, fun, axis = NA_real_, keepdims = FALSE) {
  if(is.na(axis))
    result <- fun(x)
  else {
    if(is.vector(x))
      x <- matrix(x, ncol = 1)
    result <- apply(x, axis, fun)
  }

  if(keepdims) {
    new_dim <- dim(x)
    if(is.null(new_dim))
      return(result)
    collapse <- setdiff(1:length(new_dim), axis)
    new_dim[collapse] <- 1
    dim(result) <- new_dim
  }
  result
}

####################################
#                                  #
# Utility functions for dimensions #
#                                  #
####################################
sum_dims <- function(dims) {
  if(length(dims) == 0)
    return(NULL)
  else if(length(dims) == 1)
    return(dims[[1]])

  dim <- dims[[1]]
  for(t in dims[2:length(dims)]) {
    # Only allow broadcasting for 0-D arrays or summation of scalars.
    # if(!(length(dim) == length(t) && all(dim == t)) && (!is.null(dim) && sum(dim != 1) != 0) && (!is.null(t) && sum(t != 1) != 0))
    # if(!identical(dim, t) && (!is.null(dim) && !all(dim == 1)) && (!is.null(t) && !all(t == 1)))
    if(!((length(dim) == length(t) && all(dim == t)) || all(dim == 1) || all(t == 1)))
      stop("Cannot broadcast dimensions")

    if(length(dim) >= length(t))
      longer <- dim
    else
      longer <- t

    if(length(dim) < length(t))
      shorter <- dim
    else
      shorter <- t

    offset <- length(longer) - length(shorter)
    if(offset == 0)
      prefix <- c()
    else
      prefix <- longer[1:offset]
    suffix <- c()

    if(length(shorter) > 0) {
      for(idx in length(shorter):1) {
        d1 <- longer[offset + idx]
        d2 <- shorter[idx]
        # if(!(length(d1) == length(d2) && all(d1 == d2)) && !(d1 == 1 || d2 == 1))
        if(d1 != d2 && !(d1 == 1 || d2 == 1))
          stop("Incompatible dimensions")
        if(d1 >= d2)
          new_d <- d1
        else
          new_d <- d2
        suffix <- c(new_d, suffix)
      }
    }
    dim <- c(prefix, suffix)
  }
  return(dim)
}

mul_dims_promote <- function(lh_dim, rh_dim) {
  if(is.null(lh_dim) || is.null(rh_dim) || length(lh_dim) == 0 || length(rh_dim) == 0)
    stop("Multiplication by scalars is not permitted")

  if(length(lh_dim) == 1)
    lh_dim <- c(1, lh_dim)
  if(length(rh_dim) == 1)
    rh_dim <- c(rh_dim, 1)

  lh_mat_dim <- lh_dim[(length(lh_dim)-1):length(lh_dim)]
  rh_mat_dim <- rh_dim[(length(rh_dim)-1):length(rh_dim)]
  if(length(lh_dim) > 2)
    lh_head <- lh_dim[1:(length(lh_dim)-2)]
  else
    lh_head <- c()
  if(length(rh_dim) > 2)
    rh_head <- rh_dim[1:(length(rh_dim)-2)]
  else
    rh_head <- c()

  # if(lh_mat_dim[2] != rh_mat_dim[1] || !(length(lh_head) == length(rh_head) && all(lh_head == rh_head)))
  if(lh_mat_dim[2] != rh_mat_dim[1] || !identical(lh_head, rh_head))
    stop("Incompatible dimensions")
  list(lh_dim, rh_dim, c(lh_head, lh_mat_dim[1], rh_mat_dim[2]))
}

mul_dims <- function(lh_dim, rh_dim) {
  lh_old <- lh_dim
  rh_old <- rh_dim

  promoted <- mul_dims_promote(lh_dim, rh_dim)
  lh_dim <- promoted[[1]]
  rh_dim <- promoted[[2]]
  dim <- promoted[[3]]

  # if(!(length(lh_dim) == length(lh_old) && all(lh_dim == lh_old)))
  if(!identical(lh_dim, lh_old)) {
    if(length(dim) <= 1)
      dim <- c()
    else
      dim <- dim[2:length(dim)]
  }
  # if(!(length(rh_dim) == length(rh_old) && all(rh_dim == rh_old)))
  if(!identical(rh_dim, rh_old)) {
    if(length(dim) <= 1)
      dim <- c()
    else
      dim <- dim[1:(length(dim)-1)]
  }
  return(dim)
}

###############################
#                             #
# Utility functions for signs #
#                             #
###############################
sum_signs <- function(exprs) {
  is_pos <- all(sapply(exprs, function(expr) { is_nonneg(expr) }))
  is_neg <- all(sapply(exprs, function(expr) { is_nonpos(expr) }))
  c(is_pos, is_neg)
}

mul_sign <- function(lh_expr, rh_expr) {
  # ZERO * ANYTHING == ZERO
  # NONNEGATIVE * NONNEGATIVE == NONNEGATIVE
  # NONPOSITIVE * NONNEGATIVE == NONPOSITIVE
  # NONPOSITIVE * NONPOSITIVE == NONNEGATIVE
  lh_nonneg <- is_nonneg(lh_expr)
  rh_nonneg <- is_nonneg(rh_expr)
  lh_nonpos <- is_nonpos(lh_expr)
  rh_nonpos <- is_nonpos(rh_expr)

  lh_zero <- lh_nonneg && lh_nonpos
  rh_zero <- rh_nonneg && rh_nonpos
  is_zero <- lh_zero || rh_zero

  is_pos <- is_zero || (lh_nonneg && rh_nonneg) || (lh_nonpos && rh_nonpos)
  is_neg <- is_zero || (lh_nonneg && rh_nonpos) || (lh_nonpos && rh_nonneg)
  c(is_pos, is_neg)
}

#####################################
#                                   #
# Utility functions for constraints #
#                                   #
#####################################
# Formats all the row/column cones for the solver.
format_axis <- function(t, X, axis) {
  # Reduce to norms of columns
  if(axis == 1)
    X <- lo.transpose(X)

  # Create matrices Tmat, Xmat such that Tmat*t + Xmat*X
  # gives the format for the elementwise cone constraints.
  cone_size <- 1 + nrow(X)
  terms <- list()

  # Make t_mat.
  mat_dim <- c(cone_size, 1)
  t_mat <- sparseMatrix(i = 1, j = 1, x = 1.0, dims = mat_dim)
  t_mat <- create_const(t_mat, mat_dim, sparse = TRUE)
  t_vec <- t
  if(is.null(dim(t)))   # t is scalar.
    t_vec <- lo.reshape(t, c(1,1))
  else   # t is 1-D.
    t_vec <- lo.reshape(t, c(1, nrow(t)))
  mul_dim <- c(cone_size, ncol(t_vec))
  terms <- c(terms, list(lo.mul_expr(t_mat, t_vec, mul_dim)))

  # Make X_mat.
  if(length(dim(X)) == 1)
    X <- lo.reshape(X, c(nrow(X), 1))
  mat_dim <- c(cone_size, nrow(X))
  val_arr <- rep(1.0, cone_size - 1)
  row_arr <- 2:cone_size
  col_arr <- 1:(cone_size - 1)
  X_mat <- sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = mat_dim)
  X_mat <- create_const(X_mat, mat_dim, sparse = TRUE)
  mul_dim <- c(cone_size, ncol(X))
  terms <- c(terms, list(lo.mul_expr(X_mat, X, mul_dim)))
  list(create_geq(lo.sum_expr(terms)))
}

# Formats all the elementwise cones for the solver.
format_elemwise <- function(vars_) {
  # Create matrices A_i such that 0 <= A_0*x_0 + ... + A_n*x_n
  # gives the format for the elementwise cone constraints.
  spacing <- length(vars_)

  # Matrix spaces out columns of the LinOp expressions
  var_dim <- vars_[[1]]$dim
  mat_dim <- c(spacing*var_dim[1], var_dim[1])

  mats <- lapply(0:(spacing-1), function(offset) { get_spacing_matrix(mat_dim, spacing, offset) })
  terms <- mapply(function(var, mat) { list(lo.mul_expr(mat, var)) }, vars_, mats)
  list(create_geq(lo.sum_expr(terms)))
}

# Returns a sparse matrix LinOp that spaces out an expression.
get_spacing_matrix <- function(dim, spacing, offset) {
  col_arr <- 1:dim[2]
  row_arr <- spacing*(col_arr - 1) + 1 + offset
  val_arr <- rep(1.0, dim[2])
  mat <- sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = dim)
  create_const(mat, dim, sparse = TRUE)
}

###################################
#                                 #
# Utility functions for gradients #
#                                 #
###################################
constant_grad <- function(expr) {
  grad <- list()
  for(var in variables(expr)) {
    rows <- prod(size(var))
    cols <- prod(size(expr))
    # Scalars -> 0
    if(rows == 1 && cols == 1)
      grad[[as.character(var@id)]] <- 0.0
    else
      grad[[as.character(var@id)]] <- sparseMatrix(i = c(), j = c(), dims = c(rows, cols))
  }
  grad
}

error_grad <- function(expr) {
  vars <- variables(expr)
  grad <- lapply(vars, function(var){ NA })
  names(grad) <- sapply(vars, function(var) { var@id })
  grad
}

###################
#                 #
# Power utilities #
#                 #
###################
gm <- function(t, x, y) {
  length <- size(t)
  SOC(t = reshape_expr(x+y, c(length, 1)),
      X = vstack(reshape_expr(x-y, c(1, length)), reshape_expr(2*t, c(1, length))),
      axis = 2)
}

# Form internal constraints for weighted geometric mean t <= x^p
gm_constrs <- function(t, x_list, p) {
  if(!is_weight(p))
    stop("p must be a valid weight vector")
  w <- dyad_completion(p)

  tree <- decompose(w)
  t_dim <- dim(t)
  d <- Rdictdefault(default = function(key) { new("Variable", dim = t_dim) })
  d[w] <- t

  if(length(x_list) < length(w))
    x_list <- c(x_list, list(t))

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
  for(item in tree$items) {
    elem <- item$key
    children <- item$value

    if(!any(elem == 1)) {
      for(key in list(elem, children[[1]], children[[2]])) {
        if(!is.element(key, d))
          d[key] <- d@default(key)   # Generate new value using default function
      }
      constraints <- c(constraints, list(gm(d[elem], d[children[[1]]], d[children[[2]]])))
    }
  }
  constraints
}

# TODO: For powers of 2 and 1/2 only. Get rid of this when gm_constrs is working in general.
# gm_constrs_spec <- function(t, x_list, p) {
#  list(gm(t, x_list[[1]], x_list[[2]]))
# }

# Return (t,1,x) power tuple: x <= t^(1/p) 1^(1-1/p)
pow_high <- function(p) {
  if(p <= 1) stop("Must have p > 1")
  p <- 1/gmp::as.bigq(p)
  if(1/p == as.integer(1/p))
    return(list(as.integer(1/p), c(p, 1-p)))
  list(1/p, c(p, 1-p))
}

# Return (x,1,t) power tuple: t <= x^p 1^(1-p)
pow_mid <- function(p) {
  if(p >= 1 || p <= 0) stop("Must have 0 < p < 1")
  p <- gmp::as.bigq(p)
  list(p, c(p, 1-p))
}

# Return (x,t,1) power tuple: 1 <= x^(p/(p-1)) t^(-1/(p-1))
pow_neg <- function(p) {
  if(p >= 0) stop("must have p < 0")
  p <- gmp::as.bigq(p)
  p <- p/(p-1)
  list(p/(p-1), c(p, 1-p))
}

limit_denominator <- function(num, max_denominator = 10^6) {
  # Adapted from the Python 2.7 fraction library: https://github.com/python/cpython/blob/2.7/Lib/fractions.py
  if(max_denominator < 1)
    stop("max_denominator should be at least 1")
  if(gmp::denominator(num) <= max_denominator)
    return(gmp::as.bigq(num))

  p0 <- 0
  q0 <- 1
  p1 <- 1
  q1 <- 0

  n <- as.double(gmp::numerator(num))
  d <- as.double(gmp::denominator(num))

  while(TRUE) {
    a <- floor(n/d)
    q2 <- q0 + a*q1
    if(q2 > max_denominator)
      break

    p0 <- p1
    q0 <- q1
    p1 <- p0 + a*p1
    q1 <- q2

    n <- d
    d <- n - a*d
  }

  k <- floor((max_denominator - q0)/q1)
  bound1 <- gmp::as.bigq(p0 + k*p1, q0 + k*q1)
  bound2 <- gmp::as.bigq(p1, q1)
  if(abs(bound2 - num) <= abs(bound1 - num))
    return(bound2)
  else
    return(bound1)
}

# Test if num is a positive integer power of 2
is_power2 <- function(num) {
  num > 0 && gmp::is.whole(log2(as.double(num)))
}

# Test if frac is a non-negative dyadic fraction or integer
is_dyad <- function(frac) {
  if(gmp::is.whole(frac) && frac >= 0)
    TRUE
  else if(gmp::is.bigq(frac) && frac >= 0 && is_power2(gmp::denominator(frac)))
    TRUE
  else
    FALSE
}

# Test if a vector is a valid dyadic weight vector
is_dyad_weight <- function(w) {
  is_weight(w) && all(sapply(w, is_dyad))
}

# Test if w is a valid weight vector, which consists of non-negative integer or fraction elements that sum to 1
is_weight <- function(w) {
  # if(is.matrix(w) || is.vector(w))
  #  w <- as.list(w)
  valid_elems <- rep(FALSE, length(w))
  for(i in 1:length(w))
    valid_elems[i] <- (w[i] >= 0) && (gmp::is.whole(w[i]) || gmp::is.bigq(w[i]))
  all(valid_elems) && all.equal(sum(as.double(w)), 1)
}

# Return a valid fractional weight tuple (and its dyadic completion) to represent the weights given by "a"
# When the input tuple contains only integers and fractions, "fracify" will try to represent the weights exactly
fracify <- function(a, max_denom = 1024, force_dyad = FALSE) {
  if(any(a < 0))
    stop("Input powers must be non-negative")

  if(!(gmp::is.whole(max_denom) && max_denom > 0))
    stop("Input denominator must be an integer")

  # TODO: Handle case where a contains mixture of BigQ, BigZ, and regular R numbers
  # if(is.matrix(a) || is.vector(a))
  #  a <- as.list(a)
  max_denom <- next_pow2(max_denom)
  total <- sum(a)

  if(force_dyad)
    w_frac <- make_frac(a, max_denom)
  else if(all(sapply(a, function(v) { gmp::is.whole(v) || gmp::is.bigq(v) }))) {
    w_frac <- rep(gmp::as.bigq(0), length(a))
    for(i in 1:length(a))
      w_frac[i] <- gmp::as.bigq(a[i], total)
    d <- max(gmp::denominator(w_frac))
    if(d > max_denom)
      stop("Can't reliably represent the weight vector. Try increasing 'max_denom' or checking the denominators of your input fractions.")
  } else {
    # fall through code
    w_frac <- rep(gmp::as.bigq(0), length(a))
    for(i in 1:length(a))
      w_frac[i] <- gmp::as.bigq(as.double(a[i])/total)
    if(as.double(sum(w_frac)) != 1)     # TODO: Do we need to limit the denominator with gmp?
      w_frac <- make_frac(a, max_denom)
  }
  list(w_frac, dyad_completion(w_frac))
}

# Approximate "a/sum(a)" with tuple of fractions with denominator exactly "denom"
make_frac <- function(a, denom) {
  a <- as.double(a)/sum(a)
  b <- sapply(a, function(v) { as.double(v * denom) })
  b <- as.matrix(as.integer(b))
  err <- b/as.double(denom) - a

  inds <- order(err)[1:(denom - sum(b))]
  b[inds] <- b[inds] + 1

  denom <- as.integer(denom)
  bout <- rep(gmp::as.bigq(0), length(b))
  for(i in 1:length(b))
    bout[i] <- gmp::as.bigq(b[i], denom)
  bout
}

# Return the dyadic completion of "w"
dyad_completion <- function(w) {
  # TODO: Need to handle whole decimal numbers here, e.g.
  # w <- c(1, 0, 0.0, gmp::as.bigq(0,1)) should give c(Bigq(1,1), Bigq(0,1), Bigq(0,1), Bigq(0,1))
  for(i in 1:length(w))
    w[i] <- gmp::as.bigq(w[i])
  d <- gmp::as.bigq(max(gmp::denominator(w)))

  # if extra_index
  p <- next_pow2(d)
  if(p == d)
    # the tuple of fractions is already dyadic
    return(w)
  else {
    # need to add dummy variable to represent as dyadic
    orig <- rep(gmp::as.bigq(0), length(w))
    for(i in 1:length(w))
      orig[i] <- gmp::as.bigq(w[i]*d, p)
    return(c(orig, gmp::as.bigq(p-d,p)))
  }
}

# Return the l_infinity norm error from approximating the vector a_orig/sum(a_orig) with the weight vector w_approx
approx_error <- function(a_orig, w_approx) {
  if(!all(a_orig >= 0))
    stop("a_orig must be a non-negative vector")
  if(!is_weight(w_approx))
    stop("w_approx must be a weight vector")
  if(length(a_orig) != length(w_approx))
    stop("a_orig and w_approx must have the same length")

  w_orig <- as.matrix(a_orig)/sum(a_orig)

  as.double(max(abs(w_orig - w_approx)))
}

# Return first power of 2 >= n
next_pow2 <- function(n) {
  if(n <= 0) return(1)

  ## NOTE: With CVXR 1.0, R.utils is no longer imported since log2 will work on both gmp integers
  ## and 64-bit integers
  # len <- nchar(R.utils::intToBin(n))
  # n2 <- bitwShiftL(1, len)
  # if(bitwShiftR(n2, 1) == n)
  #   return(n)
  # else
  #   return(n2)
  p <- log2(as.double(n))
  return(2^ceiling(p))
}

# Check that w_dyad is a valid dyadic completion of w
check_dyad <- function(w, w_dyad) {
  if(!(is_weight(w) && is_dyad_weight(w_dyad)))
    return(FALSE)

  if(length(w) == length(w_dyad) && all(w == w_dyad))
    # w is its own dyadic completion
    return(TRUE)

  if(length(w_dyad) == length(w) + 1) {
    if(length(w) == 0)
      return(TRUE)

    denom <- 1-w_dyad[length(w_dyad)]
    cmp <- rep(gmp::as.bigq(0), length(w_dyad)-1)
    for(i in 1:(length(w_dyad)-1))
      cmp[i] <- gmp::as.bigq(w_dyad[i], denom)
    return(all(w == cmp))
  } else
    return(FALSE)
}

# Split a tuple of dyadic rationals into two children so d_tup = 1/2*(child1 + child2)
split <- function(w_dyad) {
  # vector is all zeros with single 1, so can't be split and no children
  if(any(w_dyad == 1))
    return(list())

  bit <- gmp::as.bigq(1, 1)
  child1 <- rep(gmp::as.bigq(0), length(w_dyad))
  if(is.list(w_dyad)) {
    child2 <- rep(gmp::as.bigq(0), length(w_dyad))
    for(i in 1:length(w_dyad))
      child2[i] <- 2*w_dyad[[i]]
  } else
    child2 <- 2*w_dyad

  while(TRUE) {
    for(ind in 1:length(child2)) {
      val <- child2[ind]
      if(val >= bit) {
        child2[ind] <- child2[ind] - bit
        child1[ind] <- child1[ind] + bit
      }
      if(sum(child1) == 1)
        return(list(child1, child2))
    }
    bit <- bit/2
  }
  stop("Something wrong with input w_dyad")
}

# Recursively split dyadic tuples to produce a DAG
decompose <- function(w_dyad) {
  if(!is_dyad_weight(w_dyad))
    stop("input must be a dyadic weight vector")

  tree <- Rdict()
  todo <- list(as.vector(w_dyad))

  i <- 1
  len <- length(todo)
  while(i <= len) {
    t <- todo[[i]]
    if(!is.element(t, tree)) {
      tree[t] <- split(t)
      todo <- c(todo, tree[t])
      len <- length(todo)
    }
    i <- i + 1
  }
  tree
}

# String representation of objects in a vector
prettyvec <- function(t) {
  paste("(", paste(t, collapse = ", "), ")", sep = "")
}

# Print keys of dictionary with children indented underneath
prettydict <- function(d) {
  key_sort <- order(sapply(d$keys, get_max_denom), decreasing = TRUE)
  keys <- d$keys[key_sort]
  result <- ""
  for(vec in keys) {
    child_order <- order(sapply(d[vec], get_max_denom), decreasing = FALSE)
    children <- d[vec][child_order]
    result <- paste(result, prettyvec(vec), "\n", sep = "")
    for(child in children)
      result <- paste(result, " ", prettyvec(child), "\n", sep = "")
  }
  result
}

# Get the maximum denominator in a sequence of "BigQ" and "int" objects
get_max_denom <- function(tup) {
  denom <- rep(gmp::as.bigq(0), length(tup))
  for(i in 1:length(tup)) {
    denom[i] <- gmp::denominator(gmp::as.bigq(tup[i]))
  }
  max(denom)
}

# Return lower bound on number of cones needed to represent tuple
lower_bound <- function(w_dyad) {
  if(!is_dyad_weight(w_dyad))
    stop("w_dyad must be a dyadic weight")
  md <- get_max_denom(w_dyad)

  if(is(md, "bigq") || is(md, "bigz")) {
    ## md_int <- bit64::as.integer64(gmp::asNumeric(md))
    ## bstr <- sub("^[0]+", "", bit64::as.bitstring(md_int))   # Get rid of leading zeros
    ## lb1 <- nchar(bstr)
      ## # TODO: Should use formatBin in mpfr, but running into problems with precision
      lb1  <- log2(bit64::as.integer64(gmp::asNumeric(md)))
  } else {
      ##lb1 <- nchar(R.utils::intToBin(md))-1
      lb1  <- log2(md)
  }
  # don't include zero entries
  lb2 <- sum(w_dyad != 0) - 1
  max(lb1, lb2)
}

# Return number of cones in tree beyond known lower bounds
over_bound <- function(w_dyad, tree) {
  nonzeros <- sum(w_dyad != 0)
  return(length(tree) - lower_bound(w_dyad) - nonzeros)
}

#################
#               #
# Key utilities #
#               #
#################
Key <- function(row, col) {
    if(missing(row)) row <- NULL   # Missing row/col index implies that we select all rows/cols
    if(missing(col)) col <- NULL
    list(row = row, col = col, class = "key")
}

ku_validate_key <- function(key, dim) {   # TODO: This may need to be reassessed for consistency in handling keys.
  if(length(key) > 3)
    stop("Invalid index/slice")

  nrow <- dim[1]
  ncol <- dim[2]
  row <- ku_format_slice(key$row, nrow)
  col <- ku_format_slice(key$col, ncol)

  if(!is.null(row) && !is.null(col))
    key <- Key(row = row, col = col)
  # Change single indices for vectors into double indices
  else if(is.null(row) && !is.null(col))
    key <- Key(row = seq_len(nrow), col = col)
  else if(!is.null(row) && is.null(col))
    key <- Key(row = row, col = seq_len(ncol))
  else
    stop("A key cannot be empty")
  return(key)
}

ku_format_slice <- function(key_val, dim) {
  if(is.null(key_val))
    return(NULL)
  orig_key_val <- as.integer(key_val)

  # Return if all zero indices.
  if(all(orig_key_val == 0))
    return(orig_key_val)

  # Convert negative indices to positive indices.
  if(all(orig_key_val >= 0))
    key_val <- orig_key_val
  else if(all(orig_key_val <= 0))
    key_val <- setdiff(seq_len(dim), -orig_key_val)
  else
    stop("Only 0's may be mixed with negative subscripts")

  if(all(key_val >= 0 & key_val <= dim))
    return(key_val)
  else
    stop("Index is out of bounds for axis with size ", dim)
}

ku_slice_mat <- function(mat, key) {
  if(is.vector(mat))
    mat <- matrix(mat, ncol = 1)

  if(is.matrix(key$row) && is.null(key$col))
    select_mat  <- matrix(mat[key$row], ncol = 1)
  else if(is.null(key$row) && is.null(key$col))
    select_mat <- mat
  else if(is.null(key$row) && !is.null(key$col))
    select_mat <- mat[, key$col, drop = FALSE]
  else if(!is.null(key$row) && is.null(key$col))
    select_mat <- mat[key$row, , drop = FALSE]
  else
    select_mat <- mat[key$row, key$col, drop = FALSE]
  select_mat
}

ku_dim <- function(key, dim) {
  dims <- c()

  for(i in 1:2) {
    idx <- key[[i]]
    if(is.null(idx))
      size <- dim[i]
    else {
      selection <- (1:dim[i])[idx]
      size <- length(selection)
    }
    dims <- c(dims, size)
  }
  dims
}

