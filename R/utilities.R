#'
#' The Canonical class.
#'
#' This virtual class represents a canonical expression.
#'
#' @rdname Canonical-class
setClass("Canonical", contains = "VIRTUAL")

#' @describeIn Canonical The graph implementation of the input.
#' @return A list of \code{list(affine expression, list(constraints))}.
setMethod("canonicalize", "Canonical", function(object) { stop("Unimplemented") })

#' @describeIn Canonical The canonical form of the input.
setMethod("canonical_form", "Canonical", function(object) { canonicalize(object) })

#' @rdname expression-parts
setMethod("variables", "Canonical", function(object) { stop("Unimplemented") })

#' @rdname expression-parts
setMethod("parameters", "Canonical", function(object) { stop("Unimplemented") })

#' @rdname expression-parts
setMethod("constants", "Canonical", function(object) { stop("Unimplemented") })

#' @describeIn Canonical Information needed to reconstruct the expression aside from its arguments.
setMethod("get_data", "Canonical", function(object) { list() })

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
SOLVERS_NAME <- c(ECOS_NAME, ECOS_BB_NAME, SCS_NAME)   # TODO: Add more when we implement other solvers

# Parallel (meta) solver
PARALLEL = "parallel"

# Robust CVXOPT LDL KKT solver
ROBUST_KKTSOLVER = "robust"

# Map of constraint types
EQ_MAP = "1"
LEQ_MAP = "2"
SOC_MAP = "3"
SOC_EW_MAP = "4"
SDP_MAP = "5"
EXP_MAP = "6"
BOOL_MAP = "7"
INT_MAP = "8"

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
C_KEY = "c"
OFFSET = "offset"
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
POSITIVE = "POSITIVE"
NEGATIVE = "NEGATIVE"
UNKNOWN = "UNKNOWN"

SIGN_STRINGS = c(ZERO, POSITIVE, NEGATIVE, UNKNOWN)

################################
#                              #
# Utility functions for shapes #
#                              #
################################
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

###############################
#                             #
# Utility functions for signs #
#                             #
###############################
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

#####################################
#                                   #
# Utility functions for constraints #
#                                   #
#####################################
format_axis <- function(t, X, axis) {
  # Reduce to norms of columns
  if(axis == 1)
    X <- lo.transpose(X)
  
  # Create matrices Tmat, Xmat such that Tmat*t + Xmat*X
  # gives the format for the elementwise cone constraints.
  cone_size <- 1 + size(X)[1]
  terms <- list()
  
  # Make t_mat
  mat_size <- c(cone_size, 1)
  prod_size <- c(cone_size, size(t)[1])
  t_mat <- sparseMatrix(i = 1, j = 1, x = 1.0, dims = mat_size)
  t_mat <- create_const(t_mat, mat_size, sparse = TRUE)
  terms <- c(terms, list(lo.mul_expr(t_mat, lo.transpose(t), prod_size)))

  # Make X_mat
  mat_size <- c(cone_size, size(X)[1])
  prod_size <- c(cone_size, size(X)[2])
  val_arr <- rep(1.0, cone_size - 1)
  row_arr <- 2:cone_size
  col_arr <- 1:(cone_size - 1)
  X_mat <- sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = mat_size)
  X_mat <- create_const(X_mat, mat_size, sparse = TRUE)
  terms <- c(terms, list(lo.mul_expr(X_mat, X, prod_size)))
  list(create_geq(lo.sum_expr(terms)))
}

format_elemwise <- function(vars_) {
  # Create matrices A_i such that 0 <= A_0*x_0 + ... + A_n*x_n
  # gives the format for the elementwise cone constraints.
  spacing <- length(vars_)
  prod_size <- c(spacing*vars_[[1]]$size[1], vars_[[1]]$size[2])
  
  # Matrix spaces out columns of the LinOp expressions
  mat_size <- c(spacing*vars_[[1]]$size[1], vars_[[1]]$size[1])
  
  mats <- lapply(0:(spacing-1), function(offset) { get_spacing_matrix(mat_size, spacing, offset) })
  terms <- mapply(function(var, mat) { list(lo.mul_expr(mat, var, prod_size)) }, vars_, mats)
  list(create_geq(lo.sum_expr(terms)))
}

get_spacing_matrix <- function(size, spacing, offset) {
  require(Matrix)
  col_arr <- 1:size[2]
  row_arr <- spacing*(col_arr - 1) + 1 + offset
  val_arr <- rep(1.0, size[2])
  mat <- sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = size)
  create_const(mat, size, sparse = TRUE)
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

flatten_list <- function(x) {
  y <- list()
  rapply(x, function(x) y <<- c(y,x))
  y
}

#
# The QuadCoeffExtractor class
# 
# Given a quadratic expression of size m*n, this class extracts the coefficients 
# (Ps, Q, R) such that the (i,j) entry of the expression is given by 
# t(X) %*% Ps[[k]] %*% x + Q[k,] %*% x + R[k]
# where k = i + j*m. x is the vectorized variables indexed by id_map
# 
setClass("QuadCoeffExtractor", representation(id_map = "list", N = "numeric"))

get_coeffs.QuadCoeffExtractor <- function(object, expr) {
  if(is_constant(expr))
    return(.coeffs_constant(object, expr))
  else if(is_affine(expr))
    return(.coeffs_affine(object, expr))
  else if(is(expr, "AffineProd"))
    return(.coeffs_affine_prod(object, expr))
  else if(is(expr, "QuadOverLin"))
    return(.coeffs_quad_over_lin(object, expr))
  else if(is(expr, "Power"))
    return(.coeffs_power(object, expr))
  else if(is(expr, "MatrixFrac"))
    return(.coeffs_matrix_frac(object, expr))
  else if(is(expr, "AffAtom"))
    return(.coeffs_affine_atom(object, expr))
  else
    stop("Unknown expression type: ", class(expr))
}

.coeffs_constant.QuadCoeffExtractor <- function(object, expr) {
  if(is_scalar(expr)) {
    sz <- 1
    R <- matrix(value(expr))
  } else {
    sz <- prod(size(expr))
    R <- matrix(value(expr), nrow = sz)    # TODO: Check if this should be transposed
  }
  Ps <- lapply(1:sz, function(i) { sparseMatrix(i = c(), j = c(), dims = c(object@N, object@N)) })
  Q <- sparseMatrix(i = c(), j = c(), dims = c(sz, object@N))
  list(Ps, Q, R)
}

.coeffs_affine.QuadCoeffExtractor <- function(object, expr) {
  sz <- prod(size(expr))
  canon <- canonical_form(expr)
  prob_mat <- get_problem_matrix(list(create_eq(canon[[1]])), object@id_map)
  
  V <- prob_mat[[1]]
  I <- prob_mat[[2]]
  J <- prob_mat[[3]]
  R <- prob_mat[[4]]
  
  Q <- sparseMatrix(i = I, j = J, x = V, dims = c(sz, object@N))
  Ps <- lapply(1:sz, function(i) { sparseMatrix(i = c(), j = c(), dims = c(object@N, object@N)) })
  list(Ps, Q, as.vector(R))   # TODO: Check if R is flattened correctly
}

.coeffs_affine_prod.QuadCoeffExtractor <- function(object, expr) {
  XPQR <- .coeffs_affine(expr@args[[1]])
  YPQR <- .coeffs_affine(expr@args[[2]])
  Xsize <- size(expr@args[[1]])
  Ysize <- size(expr@args[[2]])
  
  XQ <- XPQR[[2]]
  XR <- XPQR[[3]]
  YQ <- YPQR[[2]]
  YR <- YPQR[[3]]
  m <- Xsize[1]
  p <- Xsize[2]
  n <- Ysize[2]
  
  Ps  <- list()
  Q <- sparseMatrix(i = c(), j = c(), dims = c(m*n, object@N))
  R <- rep(0, m*n)
  
  ind <- 0
  for(j in 1:n) {
    for(i in 1:m) {
      M <- sparseMatrix(i = c(), j = c(), dims = c(object@N, object@N))
      for(k in 1:p) {
        Xind <- k*m + i
        Yind <- j*p + k
        
        a <- XQ[Xind,]
        b <- XR[Xind]
        c <- YQ[Yind,]
        d <- YR[Yind]
        
        M <- M + t(a) %*% c
        Q[ind,] <- Q[ind,] + b*c + d*a
        R[ind] <- R[ind] + b*d
      }
      Ps <- c(Ps, Matrix(M, sparse = TRUE))
      ind <- ind + 1
    }
  }
  list(Ps, Matrix(Q, sparse = TRUE), R)
}

.coeffs_quad_over_lin.QuadCoeffExtractor <- function(object, expr) {
  coeffs <- .coeffs_affine(object, expr@args[[1]])
  A <- coeffs[[2]]
  b <- coeffs[[3]]
  
  P <- t(A) %*% A
  q <- Matrix(2*t(b) %*% A)
  r <- sum(b*b)
  y <- value(expr@args[[2]])
  list(list(P/y), q/y, r/y)
}

.coeffs_power.QuadCoeffExtractor <- function(object, expr) {
  if(expr@p == 1)
    return(get_coeffs(object, expr@args[[1]]))
  else if(expr@p == 2) {
    coeffs <- .coeffs_affine(object, expr@args[[1]])
    A <- coeffs[[2]]
    b <- coeffs[[3]]
    
    Ps <- lapply(1:nrow(A), function(i) { Matrix(t(A[i,]) %*% A[i,], sparse = TRUE) })
    Q <- 2*Matrix(diag(b) %*% A, sparse = TRUE)
    R <- b^2
    return(list(Ps, Q, R))
  } else
    stop("Error while processing Power")
}

.coeffs_matrix_frac <- function(object, expr) {
  coeffs <- .coeffs_affine(expr@args[[1]])
  A <- coeffs[[2]]
  b <- coeffs[[3]]
  arg_size <- size(expr@args[[1]])
  m <- arg_size[1]
  n <- arg_size[2]
  Pinv <- solve(value(expr@args[[2]]))
  
  M <- sparseMatrix(i = c(), j = c(), dims = c(object@N, object@N))
  Q <- sparseMatrix(i = c(), j = c(), dims = c(1, object@N))
  R <- 0
  
  for(i in seq(1, m*n, m)) {
    A2 <- A[i:(i+m),]
    b2 <- b[i:(i+m)]
    
    M <- M + t(A2) %*% Pinv %*% A2
    Q <- Q + 2*t(A2) %*% (Pinv %*% b2)
    R <- R + sum(b2 * (Pinv %*% b2))
  }
  list(list(Matrix(M, sparse = TRUE)), Matrix(Q, sparse = TRUE), R)
}

.coeffs_affine_atom <- function(object, expr) {
  sz <- prod(size(expr))
  Ps <- lapply(1:sz, function(i) { sparseMatrix(i = c(), j = c(), dims = c(object@N, object@N)) })
  Q <- sparseMatrix(i = c(), j = c(), dims = c(sz, object@N))
  Parg <- NA
  Qarg <- NA
  Rarg <- NA
  
  fake_args <- list()
  offsets <- list()
  offset <- 0
  for(idx in 1:length(expr@args)) {
    arg <- expr@args[[idx]]
    if(is_constant(arg))
      fake_args <- c(fake_args, list(create_const(value(arg), size(arg))))
    else {
      coeffs <- get_coeffs(object, arg)
      if(is.na(Parg)) {
        Parg <- coeffs[[1]]
        Qarg <- coeffs[[2]]
        Rarg <- coeffs[[3]]
      } else {
        Parg <- Parg + coeffs[[1]]
        Qarg <- rbind(Qarg, q)
        Rarg <- c(Rarg, r)
      }
      fake_args <- c(fake_args, list(create_var(size(arg), idx)))
      offsets[idx] <- offset
      offset <- offset + prod(size(arg))
    }
  }
  graph <- graph_implementation(expr, fake_args, size(expr), get_data(expr))
  
  # Get the matrix representation of the function
  prob_mat <- get_problem_matrix(list(create_eq(graph[[1]])), offsets)
  V <- prob_mat[[1]]
  I <- prob_mat[[2]]
  J <- prob_mat[[3]]
  R <- as.vector(prob_mat[[4]])   # TODO: Check matrix is flattened correctly
  
  # Return AX + b
  for(idx in 1:length(V)) {
    v <- V[idx]
    i <- I[idx]
    j <- J[idx]
    
    Ps[[i]] <- Ps[[i]] + v*Parg[j]
    Q[i,] <- Q[i,] + v*Qarg[j,]
    R[i] <- R[i] + v*Rarg[j]
  }
  Ps <- lapply(Ps, function(P) { Matrix(P, sparse = TRUE) })
  list(Ps, Matrix(Q, sparse = TRUE), R)
}

#
# R dictionary (inefficient, hacked together implementation).
# Note: This allows arbitrary types as both keys and values.
#
setClass("Rdict", representation(keys = "list", values = "list"), prototype(keys = list(), values = list()),
         validity = function(object) {
           if(length(object@keys) != length(object@values))
             return("Number of keys must match number of values")
           if(!all(unique(object@keys) != object@keys))
             return("Keys must be unique")
           return(TRUE)
         })

Rdict <- function(keys = list(), values = list()) {
  new("Rdict", keys = keys, values = values)
}

setMethod("$", signature(x = "Rdict"), function(x, name) {
  if(name == "items") {
    items <- rep(list(list()), length(x))
    for(i in 1:length(x)) {
      tmp <- list(key = x@keys[[i]], value = x@values[[i]])
      items[[i]] <- tmp
    }
    return(items)
  } else
    slot(x, name)
})

setMethod("length", signature(x = "Rdict"), function(x) { length(x@keys) })
setMethod("is.element", signature(el = "ANY", set = "Rdict"), function(el, set) {
  for(k in set@keys) {
    if(identical(k, el))
      return(TRUE)
  }
  return(FALSE)
})

setMethod("[", signature(x = "Rdict"), function(x, i, j, ..., drop = TRUE) {
  for(k in 1:length(x@keys)) {
    if(length(x@keys[[k]]) == length(i) && all(x@keys[[k]] == i))
      return(x@values[[k]])
  }
  stop("key ", i, " was not found")
})

setMethod("[<-", signature(x = "Rdict"), function(x, i, j, ..., value) {
  if(is.element(i, x))
    x@values[[i]] <- value
  else {
    x@keys <- c(x@keys, list(i))
    x@values <- c(x@values, list(value))
  }
  return(x)
})

setClass("Rdictdefault", representation(default = "function"), contains = "Rdict")

Rdictdefault <- function(keys = list(), values = list(), default) {
  new("Rdictdefault", keys = keys, values = values, default = default)
}

setMethod("[", signature(x = "Rdictdefault"), function(x, i, j, ..., drop = TRUE) {
  if(length(x@keys) > 0) {
    for(k in 1:length(x@keys)) {
      if(length(x@keys[[k]]) == length(i) && all(x@keys[[k]] == i))
        return(x@values[[k]])
    }
  }
  
  # TODO: Can't update in place. If key doesn't exist, want to create it with default function value.
  x@keys <- c(x@keys, list(i))
  x@values <- c(x@values, list(x@default(i)))
  return(x@values[[length(x@values)]])
})

###################
#                 #
# Power utilities #
#                 #
###################
gm <- function(t, x, y) {
  two <- create_const(2, c(1,1))
  
  length <- prod(size(t))
  SOCAxis(lo.reshape(lo.sum_expr(list(x, y)), c(length, 1)),
          lo.vstack(list(
              lo.reshape(lo.sub_expr(x, y), c(1, length)),
              lo.reshape(lo.mul_expr(two, t, size(t)), c(1, length))
            ), c(2, length)), 2)
}

# Form internal constraints for weighted geometric mean t <= x^p
gm_constrs <- function(t, x_list, p) {
  if(!is_weight(p)) stop("p must be a valid weight vector")
  w <- dyad_completion(p)

  tree <- decompose(w)
  t_size <- size(t)
  d <- Rdictdefault(default = function(key) { create_var(t_size) })
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

limit_denominator <- function(num, max_denominator = 10^6) {
  # Adapted from the Python 2.7 fraction library: https://github.com/python/cpython/blob/2.7/Lib/fractions.py
  if(max_denominator < 1)
    stop("max_denominator should be at least 1")
  if(denominator(num) <= max_denominator)
    return(as.bigq(num))
  
  p0 <- 0
  q0 <- 1
  p1 <- 1
  q1 <- 0
  
  n <- as.double(numerator(num))
  d <- as.double(denominator(num))
  
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
  bound1 <- as.bigq(p0 + k*p1, q0 + k*q1)
  bound2 <- as.bigq(p1, q1)
  if(abs(bound2 - num) <= abs(bound1 - num))
    return(bound2)
  else
    return(bound1)
}

# Test if num is a positive integer power of 2
is_power2 <- function(num) {
  num > 0 && is.whole(log2(as.double(num)))
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
  valid_elems <- rep(FALSE, length(w))
  for(i in 1:length(w))
    valid_elems[i] <- (w[i] >= 0) && (is.whole(w[i]) || is.bigq(w[i]))
  all(valid_elems) && all.equal(sum(as.double(w)), 1)
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
  
  # require(R.utils)
  # len <- nchar(intToBin(n))
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
  if(any(w_dyad == 1))
    return(list())
  
  bit <- as.bigq(1, 1)
  child1 <- rep(as.bigq(0), length(w_dyad))
  if(is.list(w_dyad)) {
    child2 <- rep(as.bigq(0), length(w_dyad))
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
  
  if(is(md, "bigq") || is(md, "bigz")) {
    require(bit64)
    md_int <- as.integer64(asNumeric(md))
    bstr <- sub("^[0]+", "", as.bitstring(md_int))   # Get rid of leading zeros
    lb1 <- nchar(bstr)
    # TODO: Should use formatBin in mpfr, but running into problems with precision
  } else {
    require(R.utils)
    lb1 <- nchar(intToBin(md))-1
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
  if(missing(row)) row <- "all"   # TODO: Missing row/col index implies that we select all rows/cols
  if(missing(col)) col <- "all"
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
