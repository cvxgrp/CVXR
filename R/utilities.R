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
SOLVERS_NAME <- c(ECOS_NAME, SCS_NAME)   # TODO: Add more when we implement other solvers

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
format_axis <- function(t, X, axis) {
  # Reduce to norms of columns
  if(axis == 1)
    X <- transpose(X)
  
  # Create matcies Tmat, Xmat such that Tmat*t + Xmat*X
  # gives the format for the elementwise cone constraints.
  cone_size <- 1 + size(X)[1]
  terms <- list()
  
  # Make t_mat
  mat_size <- c(cone_size, 1)
  prod_size <- c(cone_size, size(t)[1])
  t_mat <- sparseMatrix(i = 1, j = 1, x = 1.0, dims = mat_size)
  t_mat <- create_const(t_mat, mat_size, sparse = TRUE)
  terms <- c(terms, mul_expr(t_mat, transpose(t), prod_size))

  # Make X_mat
  mat_size <- c(cone_size, size(X)[1])
  prod_size <- c(cone_size, size(X)[2])
  val_arr <- rep(1.0, cone_size - 1)
  row_arr <- 2:cone_size
  col_arr <- 1:(cone_size-1)   # TODO: Check row_arr and col_arr indices are correct
  X_mat <- sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = mat_size)
  X_mat <- create_const(X_mat, mat_size, sparse = TRUE)
  terms <- c(terms, mul_expr(X_mat, X, prod_size))
  list(create_geq(sum_expr(terms)))
}

format_elemwise <- function(vars_) {
  # Create matcies A_i such that 0 <= A_0*x_0 + ... + A_n*x_n
  # gives the format for the elementwise cone constraints.
  spacing <- length(vars_)
  prod_size <- c(spacing * vars_[[1]]$size[1], vars_[[1]]$size[2])
  
  # Matrix spaces out columns of the LinOp expressions
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
#' Utility functions for gradients
#'
constant_grad <- function(expr) {
  grad <- list()
  for(var in variables(expr)) {
    rows <- prod(size(var))
    cols <- prod(size(expr))
    # Scalars -> 0
    if(rows == 1 && cols == 1)
      grad[var@id] = 0.0
    else
      grad[var@id] = sparseMatrix(i = c(), j = c(), dims = c(rows, cols))
  }
  grad
}

error_grad <- function(expr) {
  vars <- variables(expr)
  grad <- lapply(vars, function(var){ NA })
  names(grad) <- sapply(vars, function(var) { var@id })
  grad
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

# Form internal constraints for weighted geometric mean t <= x^p
gm_constrs <- function(t, x_list, p) {
  if(!is_weight(p)) stop("p must be a valid weight vector")
  w <- dyad_completion(p)

  # Note: This function requires a Python-style defaultdict to build and traverse a tree.
  # Currently, the closest solution is the R dict library, which is available on Github: https://github.com/mkuhn/dict
  if(!require(devtools)) install.packages("devtools")
  devtools::install_github("mkuhn/dict")
  
  tree <- decompose(w)
  # TODO: R dict must allow lists of numbers/bigq/bigz objects as keys
  d <- dict(init_keys = w, init_values = t)
  
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
  for(item in tree$items()) {
    elem <- item$key
    children <- item$value
    if(!(1 %in% elem)) {
      # If key doesn't exist, then create key and assign default value = create_var(size(key))
      if(!(elem %in% d))
        d[elem] <- create_var(size(elem))
      if(!(children[1] %in% d))
        d[children[1]] <- create_var(size(children[1]))
      if(!(children[2] %in% d))
        d[children[2]] <- create_var(size(children[2]))
      constraints <- c(constraints, list(gm(d[elem], d[children[1]], d[children[2]])))
    }
  }
  constraints
}

# Return (t,1,x) power tuple: x <= t^(1/p) 1^(1-1/p)
pow_high <- function(p) {
  if(p <= 1) stop("Must have p > 1")
  p <- 1/as.bigq(p)
  if(1/p == as.integer(1/p))
    return(list(as.integer(1/p), c(p, 1-p)))
  list(1/p, c(p, 1-p))
}

# Return (x,1,t) power tuple: t <= x^p 1^(1-p)
pow_mid <- function(p) {
  if(p >= 1 || p <= 0) stop("Must have 0 < p < 1")
  p <- as.bigq(p)
  list(p, c(p, 1-p))
}

# Return (x,t,1) power tuple: 1 <= x^(p/(p-1)) t^(-1/(p-1))
pow_neg <- function(p) {
  if(p >= 0) stop("must have p < 0")
  p <- as.bigq(p)
  p <- p/(p-1)
  list(p/(p-1), c(p, 1-p))
}

# Test if num is a positive integer power of 2
is_power2 <- function(num) {
  num > 0 && is.whole(log2(num))
}

# Test if frac is a non-negative dyadic fraction or integer
is_dyad <- function(frac) {
  if(is.whole(frac) && frac >= 0)
    TRUE
  else if(is.bigq(frac) && frac >= 0 && is_power2(denominator(frac)))
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
  valid_elems <- all(sapply(w, function(v) {
    v >= 0 && (is.whole(v) || is.bigq(v))
  }))
  valid_elems && sum(w) == 1
}

# Return a valid fractional weight tuple (and its dyadic completion) to represent the weights given by "a"
# When the input tuple contains only integers and fractions, "fracify" will try to represent the weights exactly
fracify <- function(a, max_denom = 1024, force_dyad = FALSE) {
  if(any(a < 0))
    stop("Input powers must be non-negative")
  
  if(!(is.whole(max_denom) && max_denom > 0))
    stop("Input denominator must be an integer")
  
  # TODO: Handle case where a contains mixture of BigQ, BigZ, and regular R numbers
  # if(is.matrix(a) || is.vector(a))
  #  a <- as.list(a)
  max_denom <- next_pow2(max_denom)
  total <- sum(a)
  
  if(force_dyad)
    w_frac <- make_frac(a, max_denom)
  else if(all(sapply(a, function(v) { is.whole(v) || is.bigq(v) }))) {
    w_frac <- rep(as.bigq(0), length(a))
    for(i in 1:length(a))
      w_frac[i] <- as.bigq(a[i], total)
    d <- max(denominator(w_frac))
    if(d > max_denom)
      stop("Can't reliably represent the weight vector. Try increasing 'max_denom' or checking the denominators of your input fractions.")
  } else {
    # fall through code
    w_frac <- rep(as.bigq(0), length(a))
    for(i in 1:length(a))
      w_frac[i] <- as.bigq(as.double(a[i])/total)
    if(sum(w_frac) != 1)
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
  bout <- rep(as.bigq(0), length(b))
  for(i in 1:length(b))
    bout[i] <- as.bigq(b[i], denom)
  bout
}

# Return the dyadic completion of "w"
dyad_completion <- function(w) {
  # TODO: Need to handle whole decimal numbers here, e.g.
  # w <- c(1, 0, 0.0, as.bigq(0,1)) should give c(Bigq(1,1), Bigq(0,1), Bigq(0,1), Bigq(0,1))
  for(i in 1:length(w))
    w[i] <- as.bigq(w[i])
  d <- as.bigq(max(denominator(w)))
  
  # if extra_index
  p <- next_pow2(d)
  if(p == d)
    # the tuple of fractions is already dyadic
    return(w)
  else {
    # need to add dummy variable to represent as dyadic
    orig <- rep(as.bigq(0), length(w))
    for(i in 1:length(w))
      orig[i] <- as.bigq(w[i]*d, p)
    return(c(orig, as.bigq(p-d,p)))
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
  
  require(R.utils)
  len <- nchar(intToBin(n))
  n2 <- bitwShiftL(1, len)
  if(bitwShiftR(n2, 1) == n)
    return(n)
  else
    return(n2)
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
    cmp <- rep(as.bigq(0), length(w_dyad)-1)
    for(i in 1:(length(w_dyad)-1))
      cmp[i] <- as.bigq(w_dyad[i], denom)
    return(all(w == cmp))
  } else
    return(FALSE)
}

# Split a tuple of dyadic rationals into two children so d_tup = 1/2*(child1 + child2)
split <- function(w_dyad) {
  # vector is all zeros with single 1, so can't be split and no children
  if(1 %in% w_dyad)
    return()
  
  bit <- as.bigq(1, 1)
  child1 <- rep(as.bigq(0), length(w_dyad))
  child2 <- sapply(as.list(w_dyad), function(f) { 2*f })
  
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
  
  tree <- dict()
  todo <- list(as.vector(w_dyad))
  for(t in todo) {
    if(!(t %in% tree)) {
      tree[t] <- split(t)
      todo <- c(todo, list(tree[t]))
    }
  }
  tree
}

# String representation of objects in a vector
prettyvec <- function(t) {
  paste("(", paste(t, collapse = ", "), ")", sep = "")
}

# Print keys of dictionary with children indented underneath
prettydict <- function(d) {
  # TODO: keys = sorted(list(d.keys()), key = get_max_denom, reverse = True)
  result <- ""
  for(vec in keys) {
    # TODO: children = sorted(d[vec], key = get_max_denom, reverse = False)
    result <- paste(result, prettyvec(vec), "\n", sep = "")
    for(child in children)
      result <- paste(result, " ", prettyvec(child), "\n", sep = "")
  }
  result
}

# Get the maximum denominator in a sequence of "BigQ" and "int" objects
get_max_denom <- function(tup) {
  denom <- rep(as.bigq(0), length(tup))
  for(i in 1:length(tup)) {
    denom[i] <- denominator(as.bigq(tup[i]))
  }
  max(denom)
}

# Return lower bound on number of cones needed to represent tuple
lower_bound <- function(w_dyad) {
  require(R.utils)
  if(!is_dyad_weight(w_dyad))
    stop("w_dyad must be a dyadic weight")
  md <- get_max_denom(w_dyad)
  
  lb1 <- nchar(intToBin(md))-1
  
  # don't include zero entries
  lb2 <- sum(w_dyad != 0) - 1
  max(lb1, lb2)
}

# Return number of cones in tree beyond known lower bounds
over_bound <- function(w_dyad, tree) {
  nonzeros <- sum(w_dyad != 0)
  return(length(tree) - lower_bound(w_dyad) - nonzeros)
}

#'
#' Key utilities
#'
Key <- function(row, col) {
  if(missing(row)) row <- "all"   # TODO: Missing row/col index implies that we select all rows/cols
  if(missing(col)) row <- "all"
  list(row = row, col = col, class = "key")
}

ku_validate_key <- function(key, shape) {
  nrow <- shape[1]
  ncol <- shape[2]
  
  if(length(key) > 2)
    stop("Invalid index/slice")
  else if(length(key) == 2) {
    key <- Key(row = key$row, col = key$col)
  } else if(length(key) == 1) {
    # Change single indices for vectors into double indices
    if(nrow == 1)
      key <- Key(1, key$col)
    else if(ncol == 1)
      key <- Key(key$row, 1)
    else
      stop("Invalid index/slice")
  } else
    stop("key cannot be empty")
  return(key)
}

ku_slice_mat <- function(mat, key) {
  if(key$row == "all" && key$col == "all")
    select_mat <- mat
  else if(key$row == "all")
    select_mat <- mat[, key$col]
  else if(key$col == "all")
    select_mat <- mat[key$row, ]
  else
    select_mat <- mat[key$row, key$col]
  select_mat
}

ku_size <- function(key, shape) {
  dims <- c()
  
  for(i in 1:2) {
    idx <- key[[i]]
    if(idx == "all")
      size <- shape[i]
    else {
      selection <- (1:shape[i])[idx]
      size <- length(selection)
    }
    dims <- c(dims, size)
  }
  dims
}
