##########################
#                        #
# CVXR Package Constants #
#                        #
##########################
# Constants for operators.
PLUS <- "+"
MINUS <- "-"
MUL <- "*"

# Prefix for default named variables.
VAR_PREFIX <- "var"
# Prefix for default named parameters.
PARAM_PREFIX <- "param"

# Used to overload ==.
NP_EQUAL_STR <- "equal"

# Constraint types.
EQ_CONSTR <- "=="
INEQ_CONSTR <- "<="

# Atom groups.
SOC_ATOMS <- c("GeoMean", "Pnorm", "QuadForm", "QuadOverLin", "Power")
EXP_ATOMS <- c("LogSumExp", "LogDet", "Entr", "Exp", "KLDiv", "Log", "Log1p", "Logistic")
PSD_ATOMS <- c("LambdaMax", "LambdaSumLargest", "LogDet", "MatrixFrac", "NormNuc", "SigmaMax")

# Map of constraint types
# TODO: These should be defined in a solver model.
EQ_MAP <- "1"
LEQ_MAP <- "2"
SOC_MAP <- "3"
SOC_EW_MAP <- "4"
PSD_MAP <- "5"
EXP_MAP <- "6"
BOOL_MAP <- "7"
INT_MAP <- "8"

# Keys in the dictionary of cone dimensions.
# Cone dims are now defined in matrix stuffing modules rather than the solver module.

## EQ_DIM <- "f"   # SCS 3.0 changes this to "z"
EQ_DIM <- "z"
LEQ_DIM <- "l"
SOC_DIM <- "q"
PSD_DIM <- "s"
EXP_DIM <- "ep"

# Keys for non-convex constraints.
BOOL_IDS <- "bool_ids"
BOOL_IDX <- "bool_idx"
INT_IDS <- "int_ids"
INT_IDX <- "int_idx"

# Parametrized problem.
PARAM_PROB <- "param_prob"

# Keys for problem data dict.
C_KEY <- "c"
OFFSET <- "offset"
P_KEY <- "P"
Q_KEY <- "q"
A_KEY <- "A"
B_KEY <- "b"
G_KEY <- "G"
H_KEY <- "h"
F_KEY <- "F"
DIMS <- "dims"
BOOL_IDX <- "bool_vars_idx"
INT_IDX <- "int_vars_idx"

# Keys for curvature and sign
CONSTANT <- "CONSTANT"
AFFINE <- "AFFINE"
CONVEX <- "CONVEX"
CONCAVE <- "CONCAVE"
QUASILINEAR <- "QUASILINEAR"
QUASICONVEX <- "QUASICONVEX"
QUASICONCAVE <- "QUASICONCAVE"
LOG_LOG_CONSTANT <- "LOG-LOG CONSTANT"
LOG_LOG_AFFINE <- "LOG-LOG AFFINE"
LOG_LOG_CONVEX <- "LOG-LOG CONVEX"
LOG_LOG_CONCAVE <- "LOG-LOG CONCAVE"
ZERO <- "ZERO"
NONNEG <- "NONNEGATIVE"
NONPOS <- "NONPOSITIVE"
UNKNOWN <- "UNKNOWN"

# Canonicalization backends.
RUST_CANON_BACKEND <- "RUST"
CPP_CANON_BACKEND <- "CPP"
DEFAULT_CANON_BACKEND <- CPP_CANON_BACKEND

# Numerical tolerances.
EIGVAL_TOL <- 1e-10
PSD_NSD_PROJECTION_TOL <- 1e-8
GENERAL_PROJECTION_TOL <- 1e-10
SPARSE_PROJECTION_TOL <- 1e-10
ATOM_EVAL_TOL <- 1e-4

# DPP is slow when total size of parameters exceed this threshold.
PARAM_THRESHOLD <- 1e4   # TODO: Should we reduce this?

# TODO: Can we set the number of threads to use during compilation like in CVXPY?

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
squeezed <- function(dim) {
  return(dim[dim != 1])
}

sum_dims <- function(dims) {
  # Give the dimensions resulting from summing a list of dimensions.
  if(is.null(dims) || length(dims) == 0)
    return(NULL)   # Should I return NULL or 0 for scalars?
  else if(length(dims) == 1)
    return(dims[[1]])

  dim <- dims[[1]]
  for(t in dims[2:length(dims)]) {
    # Only allow broadcasting for 0-D arrays or summation of scalars.
    if(length(dim) != t && length(squeezed(dim)) != 0 && length(squeezed(t)) != 0)
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
        if(d1 != d2 && !(d1 == 1 || d2 == 1))
          stop("Incompatible dimensions")
        if(d1 >= d2)
          new_dim <- d1
        else
          new_dim <- d2
        suffix <- c(new_dim, suffix)
      }
    }
    dim <- c(prefix, suffix)
  }
  return(dim)
}

mul_dims_promote <- function(lh_dim, rh_dim) {
  # Promotes dims as necessary and returns promoted dim of product.
  # If lh_dim is of length one, prepend a one to it.
  # If rh_dim is of length one, append a one to it.

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
  # Give the dim resulting from multiplying two dims.
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

size_from_dim <- function(dim) {
  # Compute the size of a given shape by multiplying the sizes of each axis.
  # This is a replacement for as.integer(prod(dim)), which is much slower for
  # small arrays than the implementation below.
  return(Reduce("*", dim, 1))
}

###############################
#                             #
# Utility functions for signs #
#                             #
###############################
sum_signs <- function(exprs) {
  # Give the sign resulting from summing a list of expressions.
  is_pos <- all(sapply(exprs, is_nonneg))
  is_neg <- all(sapply(exprs, is_nonpos)
  c(is_pos, is_neg)
}

mul_sign <- function(lh_expr, rh_expr) {
  # Give the sign resulting from multiplying two expressions.
  # ZERO * ANYTHING == ZERO
  # POSITIVE * POSITIVE == POSITIVE
  # NEGATIVE * POSITIVE == NEGATIVE
  # NEGATIVE * NEGATIVE == POSITIVE

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
  col_arr <- seq_len(dim[2])
  row_arr <- spacing * (col_arr - 1) + 1 + offset
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
  # Returns the gradient of constant terms in an expression.
  # Matrix expressions are vectorized, so the gradient is a matrix.
  grad <- list()
  for(var in variables(expr)) {
    rows <- size(var)
    cols <- size(expr)
    # Scalars -> 0
    if(rows == 1 && cols == 1)
      grad[[as.character(id(var))]] <- 0.0
    else
      grad[[as.character(id(var))]] <- sparseMatrix(i = c(), j = c(), dims = c(rows, cols), repr = "C")
  }
  return(grad)
}

error_grad <- function(expr) {
  # Returns a gradient of all NAs.
  vars <- variables(expr)
  grad <- as.list(rep(NA_real_, length(vars)))
  names(grad) <- sapply(vars, id)
  return(grad)
}

############################
#                          #
# Linear algebra utilities #
#                          #
############################
orth <- function(V, tol = 1e-12, ...) {
  # Return a matrix whose columns are an orthonormal basis for range(V).
  res <- qr(V, tol = tol, ...)
  QR <- res$qr   # Q in upper triangle, R in lower triangle.
  rank <- res$rank
  Q_full <- QR
  Q_full[lower.tri(QR)] <- 0
  Q <- matrix(Q_full[,1:rank], nrow = nrow(V), ncol = rank)   # Ensure 2-dimensional.
  return(Q)
}

onb_for_orthogonal_complement <- function(V) {
  # Let U = the orthogonal complement of range(V).
  # This function returns an array Q whose columns are an orthonormal basis for U.
  # It requires that dim(U) > 0.
  n <- nrow(V)
  Q1 <- orth(V)
  rank <- ncol(Q1)

  if(n <= rank)
    stop("Must have n > rank")

  if(is.complex(V))
    P <- diag(n) - Q1 %*% t(Conj(Q1))
  else
    P <- diag(n) - Q1 %*% t(Q1)
  Q2 <- orth(P)
  return(Q2)
}

is_psd_within_tol <- function(A, tol) {
  # Return TRUE if we can certify that A is PSD (up to tolerance "tol").
  #
  # First we check if A is PSD according to the Gershgorin Circle Theorem.
  #
  # If Gershgorin is inconclusive, then we use an iterative method to estimate
  # extremal eigenvalues of certain shifted versions of A. The shifts are chosen
  # so that the signs of those eigenvalues tell us the signs of the eigenvalues of A.
  #
  # If there are numerical issues then it's possible that this function returns
  # FALSE even when A is PSD. If you know that you're in that situation, then
  # you should replace A by PSDWrap(A).
  #
  # Parameters
  # -----------
  # A = Symmetric (or Hermitian) dense or sparse matrix.
  # tol = Nonnegative floating point variable. Something very small, like 1e-10.

  if(gershgorin_psd_check(A, tol))
    return(TRUE)

  # Returns the eigenvalue w[i] of A where 1/(w[i] - sigma) is minimized.
  # If A - sigma*I is PSD, then w[i] should be equal to the largest eigenvalue of A.
  # If A - sigma*I is not PSD, then w[i] should be the largest eigenvalue of A where w[i] - sigma < 0.
  # We should only call this function with sigma < 0. In this case, if A - sigma*I is not PSD,
  # then A is not PSD, and w[i] < -abs(sigma) is a negative eigenvalue of A.
  # If A - sigma*I is PSD, then we obviously have that the smallest eigenvalue of A is >= sigma.
  SA_eigsh <- function(sigma) {
    require(RSpectra)
    res <- RSpectra::eigs_sym(A, k = 1, which = "SA", sigma = sigma, opts = list(retvec = FALSE))
    return(res$values)
  }

  ev <- NA_real_
  tryCatch({
    ev <- SA_eigsh(-tol)   # Might return NA, or raise an exception.
  }, finally = {
    if(all(is.na(ev))) {
      # Will be NA if A has an eigenvalue which is exactly -tol.
      # (We might also hit this code block for other reasons).
      temp <- tol - .Machine$double.eps
      ev <- SA_eigsh(-temp)
    }
  })
  return(all(ev >= -tol))
}

gershgorin_psd_check <- function(A, tol) {
  # Use the Gershgorin Circle Theorem
  #
  # https://en.wikipedia.org/wiki/Gershgorin_circle_theorem
  #
  # As a sufficient condition for A being PSD with tolerance "tol".
  #
  # The computational complexity of this function is O(nnz(A)).
  #
  # Parameters
  # -----------
  # A = Symmetric (or Hermitian) dense or sparse matrix.
  # tol = Nonnegative floating point variable. Something very small, like 1e-10.

  if(is(A, "sparseMatrix")) {
    d <- diag(A)
    if(any(d < -tol))
      return(FALSE)
    A_shift <- A - sparseMatrix(i = 1:length(d), j = 1:length(d), x = d)
    radii <- apply(A_shift, MARGIN = 2, sum)
    return(all(diag - radii >= -tol))
  } else if(is.numeric(A)) {
    d <- diag(A)
    if(any(diag < -tol))
      return(FALSE)
    A_shift <- A - diag(d)
    A_shift <- abs(A_shift)
    radii <- apply(A_shift, MARGIN = 2, sum)
    return(all(diag - radii >= -tol))
  } else
    stop("A must be a sparse or dense numeric matrix")
}

#########################
#                       #
# Perspective utilities #
#                       #
#########################
form_cone_constraint <- function(z, constraint) {
  # Given a constraint represented as Ax + b in K for K a CVXR cone,
  # return an instantiated CVXR constraint.
  if(is(constraint, "SOC")) {
    # TODO: Figure out how to instantiate Ax + b in SOC where we know which
    # lines from our ultimate A_pers(x,t,s) + b in K times ... correspond to
    # this constraint.
    return(SOC(t = z[1], X = z[2:nrow(z)]))
  } else if(is(constraint, "NonNegConstraint"))
    return(NonNegConstraint(z))
  else if(is(constraint, "ExpCone")) {
    n <- nrow(z)
    if(!(length(dim(z)) == 1 || ncol(z) == 1))
      stop("z must be a vector or matrix with a single column")
    if(n %% 3 != 0)   # We think this is how the exponential cone works.
      stop("n needs to be a multiple of 3")
    step <- floor(n/3)
    return(ExpCone(z[1:step], z[(step+1):(n-step)], z[(n-step+1):n]))
  } else if(is(constraint, "ZeroConstraint"))
    return(ZeroConstraint(z))
  else if(is(constraint, "PSDConstraint")) {
    N <- nrow(z)
    n <- as.integer(N^0.5)
    if(N != n^2)
      stop("Argument is not a vectorized square matrix")
    z_mat <- Reshape(z, c(n, n))
    return(PSDConstraint(z_mat))   # Do we need constraint_id?
  } else if(is(constraint, "PowCone3D"))
    stop("Unimplemented")
  else
    stop("Unimplemented")
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

gm_constrs <- function(t, x_list, p) {
  # Form internal constraints for weighted geometric mean t <= x^p
  # t <= x[1]^p[1] * x[2]^p[2] * ... * x[n]^p[n],
  # where x and t can either be scalar or matrix variables.

  if(!is_weight(p))
    stop("p must be a valid weight vector")
  w <- dyad_completion(p)

  tree <- decompose(w)
  t_dim <- dim(t)
  d <- Rdictdefault(default = function(key) { new("Variable", dim = t_dim) })
  d[w] <- t

  long_w <- length(w) - length(x_list)
  if(long_w > 0)
    x_list <- c(x_list, as.list(rep(t, long_w)))

  if(length(x_list) != length(w))
    stop("Expected length of x_list to be equal to length of w, but got ", length(x_list), " != ", length(w))

  for(i in seq_along(x_list)) {
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
  return(constraints)
}

# TODO: For powers of 2 and 1/2 only. Get rid of this when gm_constrs is working in general.
# gm_constrs_spec <- function(t, x_list, p) {
#  list(gm(t, x_list[[1]], x_list[[2]]))
# }

pow_high <- function(p) {
  # Return (t,1,x) power tuple:
  # x <= t^(1/p) 1^(1-1/p).
  # User wants the epigraph variable t.

  if(p <= 1)
    stop("Must have p > 1")
  p <- 1/gmp::as.bigq(p)
  if(1/p == as.integer(1/p))
    return(list(as.integer(1/p), c(p, 1-p)))
  return(list(1/p, c(p, 1-p)))
}

pow_mid <- function(p) {
  # Return (x,1,t) power tuple:
  # t <= x^p 1^(1-p).
  # User wants the epigraph variable t.

  if(p >= 1 || p <= 0)
    stop("Must have 0 < p < 1")
  p <- gmp::as.bigq(p)
  return(list(p, c(p, 1-p)))
}

pow_neg <- function(p) {
  # Return (x,t,1) power tuple:
  # 1 <= x^(p/(p-1)) t^(-1/(p-1)).
  # User wants the epigraph variable t.

  if(p >= 0)
    stop("must have p < 0")
  p <- gmp::as.bigq(p)
  p <- p/(p-1)
  return(list(p/(p-1), c(p, 1-p)))
}

limit_denominator <- function(num, max_denominator = 10^6) {
  # Closest fraction to self with denominator at most max_denominator.
  # Adapted from the Python fractions library: https://github.com/python/cpython/blob/main/Lib/fractions.py

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

  # Determine which of the candidates (p0+k*p1)/(q0+k*q1) and p1/q1 is
  # closer to self. The distance between them is 1/(q1*(q0+k*q1)), while
  # the distance from p1/q1 to self is d/(q1*self._denominator). So we
  # need to compare 2*(q0+k*q1) with self._denominator/d.
  denom <- as.double(gmp::denominator(num))
  if((2*d*(q0 + k*q1)) <= denom)
    return(gmp::as.bigq(p1, q1))
  else
    return(gmp::as.bigq(p0 + k*p1, q0 + k*q1))
}

is_power2 <- function(num) {
  # Test if num is a positive integer power of 2.
  # Note: Unlike Python, this uses the actual value of the number, so is_power2(1.0) returns TRUE.
  return(gmp::is.whole(num) && num > 0 && !bitwAnd(num, num - 1))
}

is_dyad <- function(frac) {
  # Test if frac is a non-negative dyadic fraction or integer.
  if(gmp::is.whole(frac) && frac >= 0)
    return(TRUE)
  else if(gmp::is.bigq(frac) && frac >= 0 && is_power2(gmp::denominator(frac)))
    return(TRUE)
  else
    return(FALSE)
}

is_dyad_weight <- function(w) {
  # Test if a vector is a valid dyadic weight vector.
  # w must be nonnegative, sum to 1, and have integer or dyadic fractional elements.
  return(is_weight(w) && all(sapply(w, is_dyad)))
}

is_weight <- function(w) {
  # Test if w is a valid weight vector.
  # w must have nonnegative integer or fraction elements, and sum to 1.

  # if(is.matrix(w) || is.vector(w))
  #   w <- as.list(w)

  valid_elems <- rep(FALSE, length(w))
  for(i in seq_along(w))
    valid_elems[i] <- (w[i] >= 0) && (gmp::is.whole(w[i]) || gmp::is.bigq(w[i]))
  valid_elems <- all(valid_elems)

  # return(all(valid_elems) && all.equal(sum(as.double(w)), 1))
  return(valid_elems && sum(w) == 1)
}

fracify <- function(a, max_denom = 1024, force_dyad = FALSE) {
  # Return a valid fractional weight tuple (and its dyadic completion) to represent the weights given by "a"
  # When the input tuple contains only integers and fractions, "fracify" will try to represent the weights exactly

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
  # Check if the key is a valid index.
  if(length(key) > 3)
    stop("Invalid index/slice")

  nrow <- dim[1]
  ncol <- dim[2]
  row <- ku_format_slice(key$row, nrow)
  col <- ku_format_slice(key$col, ncol)

  # Change single indices for vectors into double indices
  if(!is.null(row) && !is.null(col))
    key <- Key(row = row, col = col)
  else if(is.null(row) && !is.null(col))
    key <- Key(row = seq_len(nrow), col = col)
  else if(!is.null(row) && is.null(col))
    key <- Key(row = row, col = seq_len(ncol))
  else
    stop("A key cannot be empty")
  return(key)
}

ku_format_slice <- function(key_val, dim, axis) {
  # Converts part of a key into a slice with a start and step.
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

###########################
#                         #
# Replace quadratic forms #
#                         #
###########################

# These functions are used to handle both S4 class and list inputs to the
# utility functions for replacing quadratic forms (necessary in e.g. coeff_quad_form).
# TODO: Change LinOp representation from list to S4 class so we can get rid of
# these extra checks.
get_obj_slot(obj, name) {
  if(is.list(obj))
    return(obj[[name]])
  else
    return(slot(obj, name))
}

set_obj_slot(obj, name, idx, val) {
  if(is.list(obj))
    obj[[name]][[idx]] <- val
  else
    slot(obj, name)[[idx]] <- val
  return(obj)
}

replace_quad_forms <- function(expr, quad_forms) {
  # for(idx in seq_along(expr@args)) {
  #   arg <- expr@args[[idx]]
  for(idx in seq_along(get_obj_slot(expr, "args"))) {
    arg <- get_obj_slot(expr, "args")[[idx]]
    if(is(arg, "SymbolicQuadForm"), || is(arg, "QuadForm")) {
      tmp <- replace_quad_form(expr, idx, quad_forms)
      expr <- tmp[[1]]
      quad_forms <- tmp[[2]]
    } else {
      tmp <- replace_quad_forms(arg, quad_forms)
      arg <- tmp[[1]]
      quad_forms <- tmp[[2]]
    }
  }
  return(list(expr, quad_forms))
}

replace_quad_form <- function(expr, idx, quad_forms) {
  # quad_form <- expr@args[[idx]]
  quad_form <- get_obj_slot(expr, "args")[[idx]]
  placeholder <- new("Variable", dim = dim(quad_form), var_id = id(quad_form))
  # expr@args[[idx]] <- placeholder
  expr <- set_obj_slot(expr, "args", idx, placeholder)
  placeholder_id_char <- as.character(id(placeholder))
  quad_forms[[placeholder_id_char]] <- list(expr, idx, quad_form)
  return(list(expr, quad_forms))
}

restore_quad_forms <- function(expr, quad_forms) {
  # TODO: Check recursion is handled correctly by returning expr. (CVXPY modifies expr@args in place).
  # for(idx in seq_along(expr@args)) {
  #   arg <- expr@args[[idx]]
  for(idx in seq_along(get_obj_slot(expr, "args"))) {
    arg <- get_obj_slot(expr, "args")
    arg_id_char <- as.character(id(arg))
    if(is(arg, "Variable") && arg_id_char %in% names(quad_forms)) {
      # expr@args[[idx]] <- quad_forms[[arg_id_char]][[3]]
      expr <- set_obj_slot(expr, "args", idx, quad_forms[[arg_id_char]][[3]])
    } else {
      # restore_quad_forms(arg, quad_forms)
      # expr@args[[idx]] <- restore_quad_forms(arg, quad_forms)
      expr <- set_obj_slot(expr, "args", idx, restore_quad_forms(arg, quad_forms))
    }
  }
  return(expr)
}

#######################
#                     #
# Debugging utilities #
#                     #
#######################
node_count <- function(expr) {
  # Return node count for the expression/constraint.
  if("args" %in% slotNames(expr)) {
    return(1 + sum(sapply(expr@args, node_count)))
  } else
    return(1)
}

build_non_disciplined_error_msg <- function(problem, discipline_type) {
  prop_name <- NA_character_
  prefix_conv <- ""
  if(discipline_type == "DCP")
    prop_name <- "is_dcp"
  else if(discipline_type == "DGP") {
    prop_name <- "is_dgp"
    prefix_conv <- "log_log_"
  } else
    stop("Unknown discipline type ", discipline_type)

  find_non_prop_leaves <- function(expr, res = NULL) {
    if(is.null(res))
      res <- c()
    # if((is.null(expr@args) || length(expr@args) == 0) && slot(expr, prop_name)())
    if((is.null(expr@args) || length(expr@args) == 0) && do.call(prop_name, expr))
      return(res)

    # if(!slot(expr, prop_name)() && all(sapply(expr@args, function(child) { slot(child, prop_name)() } ))) {
    if(!do.call(prop_name, expr) && all(sapply(expr@args, function(child) { do.call(prop_name, child) } ))) {
      str_expr <- as.character(expr)
      if(discipline_type == "DGP" && is(expr, "Variable"))
        str_expr <- paste(str_expr, "<-- needs to be declared positive")
      res <- c(res, str_expr)
    }

    for(child in expr@args)
      res <- find_non_prop_leaves(child, res)
    return(res)
  }

  # if(!slot(problem@objective, prop_name)()) {
  if(!do.call(prop_name, problem@objective)) {
    non_disciplined_leaves <- find_non_prop_leaves(problem@objective@expr)
    if(length(non_disciplined_leaves) > 0)
      msg <- paste("The objective is not ", discipline_type, ". Its following subexpressions are not:", sep = "")
    else {
      convex_str <- paste(prefix_conv, "convex", sep = "")
      concave_str <- paste(prefix_conv, "concave", sep = "")
      # fun_attr_check <- slot(problem@objective@args[[1]], paste("is_", convex_str, sep = ""))()
      fun_attr_check <- do.call(paste("is_", convex_str, sep = ""), problem@objective@args[[1]])
      msg <- paste("The objective is not ", discipline_type, ", even though each sub-expression is.\n",
                   "You are trying to ", problem@objective@NAME, " a function that is ",
                   ifelse(fun_attr_check, convex_str, concave_str), sep = "")
    }

    for(expr in non_disciplined_leaves)
      msg <- paste(msg, as.character(expr), sep = "\n")
    return(msg)
  }

  disciplined_mask <- sapply(problem@constraints, is_dcp)
  not_disciplined_constraints <- problem@constraints[!disciplined_mask]

  msg <- paste("The following constraints are not", discipline_type)
  for(expr in not_disciplined_constraints) {
    msg <- paste(msg, "\n", as.character(expr), ", because the following subexpressions are not:", sep = "")
    non_disciplined_leaves <- find_non_prop_leaves(expr)
    for(subexpr in non_disciplined_leaves)
      msg <- paste(msg, "\n|-- ", as.character(subexpr), sep = "")
  }
  return(msg)
}

