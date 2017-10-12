# Scalar functions
geo_mean <- GeoMean
harmonic_mean <- HarmonicMean
lambda_max <- LambdaMax
lambda_min <- LambdaMin
lambda_sum_largest <- LambdaSumLargest
lambda_sum_smallest <- LambdaSumSmallest
log_det <- LogDet
log_sum_exp <- LogSumExp
matrix_frac <- MatrixFrac
max_entries <- MaxEntries
min_entries <- MinEntries

#'
#' Mixed Norm
#' 
#' The \eqn{l_{p,q}} norm \eqn{l_{p,q}(x) = \left(\sum_{i=1}^n (\sum_{j=1}^m |x_{i,j}|)^{q/p}\right)^{1/q}}.
#' 
#' @param X An \S4class{Expression}, vector, or matrix.
#' @param p The type of inner norm.
#' @param q The type of outer norm.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \S4class{Expression} representing the \eqn{l_{p,q}} norm.
#' @aliases mixed_norm
#' @export
mixed_norm <- MixedNorm


#'
#' 1-Norm
#' 
#' The 1-norm \eqn{\|x\|_1 = \sum_{i=1}^n |x_i|}.
#' 
#' @param x An \S4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \S4class{Expression} representing the 1-norm.
#' @aliases norm1
#' @export
norm1 <- Norm1

#'
#' Euclidean Norm
#' 
#' The Euclidean norm \eqn{\|x\|_2 = \left(\sum_{i=1}^n x_i^2\right)^{1/2}}.
#' 
#' @param x An \S4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \S4class{Expression} representing the Euclidean norm.
#' @aliases norm2
norm2 <- Norm2

#'
#' Infinity-Norm
#' 
#' The \eqn{\infty}}-norm \eqn{\|x\|_{\infty} = \max_{i=1,\ldots,n} |x_i|}.
#' 
#' @param x An \S4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \S4class{Expression} representing the infinity-norm.
#' @aliases norm_inf
#' @export
norm_inf <- NormInf

#'
#' Nuclear Norm
#' 
#' The nuclear norm, i.e. sum of the singular values of a matrix.
#' 
#' @param A An \S4class{Expression} or matrix.
#' @return An \S4class{Expression} representing the nuclear norm.
#' @aliases norm_nuc
#' @export
norm_nuc <- NormNuc
p_norm <- Pnorm

#'
#' Quadratic Form
#'
#' The quadratic form \eqn{x^TPx}.
#'
#' @param x An \S4class{Expression} or vector.
#' @param P An \S4class{Expression} or matrix.
#' @return An \S4class{Expression} representing the quadratic form.
#' @aliases quad_form
#' @export
quad_form <- QuadForm

#'
#' Quadratic over Linear
#'
#' The sum of squared entries in X divided by a scalar y \eqn{\sum_{i,j} X_{i,j}^2/y}.
#'
#' @param x An \S4class{Expression}, vector, or matrix.
#' @param y A scalar \S4class{Expression} or numeric constant.
#' @return An \S4class{Expression} representing the quadratic over linear function value.
#' @aliases quad_over_lin
#' @export
quad_over_lin <- QuadOverLin

#'
#' Sum of Entries
#'
#' The sum of entries in a vector or matrix.
#'
#' @slot expr An \S4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code[NA}.
#' @return An \S4class{Expression} representing the sum of entries.
#' @aliases sum_entries
#' @export
sum_entries <- SumEntries

#'
#' Sum of Largest Values
#' 
#' The sum of the largest k values of a vector or matrix.
#' 
#' @param x An \S4class{Expression}, vector, or matrix.
#' @param k The number of largest values to sum over.
#' @return An \S4class{Expression} representing the sum of the largest k values.
#' @aliases sum_largest
#' @export
sum_largest <- SumLargest

#'
#' Sum of Smallest Values
#' 
#' The sum of the smallest k values of a vector or matrix.
#' 
#' @param x An \S4class{Expression}, vector, or matrix.
#' @param k The number of smallest values to sum over.
#' @return An \S4class{Expression} representing the sum of the smallest k values.
#' @aliases sum_smallest
#' @export
sum_smallest <- SumSmallest

#'
#' Sum of Squares
#' 
#' The sum of the squared entries in a vector or matrix.
#' 
#' @param expr An \S4class{Expression}, vector, or matrix.
#' @return An \S4class{Expression} representing the sum of squares.
#' @aliases sum_squares
#' @export
sum_squares <- SumSquares

#'
#' Matrix Trace
#'
#' The sum of the diagonal entries in a matrix.
#'
#' @param expr An \S4class{Expression} or matrix.
#' @return An \S4class{Expression} representing the trace.
#' @aliases matrix_trace
matrix_trace <- Trace

#'
#' Total Variation
#' 
#' The total variation of a vector, matrix, or list of matrices. Uses L1 norm of discrete gradients for vectors and L2 norm of discrete gradients for matrices.
#' 
#' @param value An \S4class{Expression}, vector, or matrix.
#' @param ... \S4class{Expression} objects or numeric constants that extend the third dimension of value.
#' @return An \S4class{Expression} representing the total variation.
#' @aliases tv
#' @export
tv <- TotalVariation

max.Expression <- function(..., na.rm = FALSE) {
  if(na.rm)
    warning("na.rm is unimplemented for Expression objects")
  
  vals <- list(...)
  is_expr <- sapply(vals, function(v) { is(v, "Expression") })
  max_args <- lapply(vals[is_expr], function(expr) { MaxEntries(expr) })
  if(!all(is_expr)) {
    max_num <- max(sapply(vals[!is_expr], function(v) { max(v, na.rm = na.rm) }))
    max_args <- c(max_args, max_num)
  }
  .MaxElemwise(args = max_args)
}

min.Expression <- function(..., na.rm = FALSE) {
  if(na.rm)
    warning("na.rm is unimplemented for Expression objects")
  
  vals <- list(...)
  is_expr <- sapply(vals, function(v) { is(v, "Expression") })
  min_args <- lapply(vals[is_expr], function(expr) { MinEntries(expr) })
  if(!all(is_expr)) {
    min_num <- min(sapply(vals[!is_expr], function(v) { min(v, na.rm = na.rm) }))
    min_args <- c(min_args, min_num)
  }
  min_args <- lapply(min_args, function(arg) { -as.Constant(arg) })
  -.MaxElemwise(args = min_args)
}

setMethod("norm", signature(x = "Expression", type = "character"), function(x, type) {
  x <- as.Constant(x)
  type <- substr(type, 1, 1)
  
  # Norms for scalars same as absolute value
  if(type %in% c("O", "o", "1"))                  # Maximum absolute column sum
    MaxEntries(Pnorm(x = x, p = 1, axis = 2))
  else if(type %in% c("I", "i"))                  # Maximum absolute row sum
    MaxEntries(Pnorm(x = x, p = 1, axis = 1))
  else if(type %in% c("E", "e", "F", "f"))        # Frobenius norm (Euclidean norm if x is treated as a vector)
    Pnorm(x = x, p = 2, axis = NA_real_)
  else if(type %in% c("M", "m"))                  # Maximum modulus (absolute value) of all elements in x
    MaxEntries(Abs(x = x))
  else if(type == "2")                            # Spectral norm (largest singular value of x)
    SigmaMax(A = x)
  else
    stop("argument type[1]='", type, "' must be one of 'M','1','O','I','F' or 'E'")
})

sum.Expression <- function(..., na.rm = FALSE) {
  if(na.rm)
    warning("na.rm is unimplemented for Expression objects")
  
  vals <- list(...)
  is_expr <- sapply(vals, function(v) { is(v, "Expression") })
  sum_expr <- lapply(vals[is_expr], function(expr) { SumEntries(expr = expr) })
  if(all(is_expr))
    Reduce("+", sum_expr)
  else {
    sum_num <- sum(sapply(vals[!is_expr], function(v) { sum(v, na.rm = na.rm) }))
    Reduce("+", sum_expr) + sum_num
  }
}

mean.Expression <- function(x, trim = 0, na.rm = FALSE, ...) {
  if(na.rm)
    stop("na.rm is unimplemented for Expression objects")
  if(trim != 0)
    stop("trim is unimplemented for Expression objects")
  SumEntries(expr = x) / prod(size(x))
}

# Elementwise functions
entr <- Entr
huber <- Huber
inv_pos <- InvPos
kl_div <- KLDiv
logistic <- Logistic
max_elemwise <- MaxElemwise
min_elemwise <- MinElemwise
mul_elemwise <- MulElemwise
neg <- Neg
pos <- Pos
power <- Power
scalene <- Scalene
square <- Square

setMethod("abs", "Expression", function(x) { Abs(x = x) })
setMethod("exp", "Expression", function(x) { Exp(x = x) })
setMethod("log", "Expression", function(x, base = exp(1)) { Log(x = x)/log(base) })
setMethod("log10", "Expression", function(x) { log(x, base = 10) })
setMethod("log2", "Expression", function(x) { log(x, base = 2) })
setMethod("log1p", "Expression", function(x) { Log1p(x = x) })
setMethod("sqrt", "Expression", function(x) { Sqrt(x = x) })

# Matrix/vector operations
affine_prod <- AffineProd
bmat <- Bmat
conv <- Conv
cum_sum <- CumSum
hstack <- HStack
kron <- Kron
reshape_expr <- Reshape

#'
#' Maximum Singular Value
#' 
#' The maximum singular value of a matrix.
#' 
#' @param A An \S4class{Expression} or matrix.
#' @return An \S4class{Expression} representing the maximum singular value.
#' @aliases sigma_max
#' @export
sigma_max <- SigmaMax

upper_tri <- UpperTri
vec <- Vec
vstack <- VStack

setMethod("cumsum", "Expression", function(x) { CumSum(expr = Vec(x)) })   # Flatten matrix in column-major order to match R's behavior
setMethod("diag", signature(x = "Expression"), function(x, nrow, ncol) {
  if(nargs() == 1L)
    Diag(x)
  else if(is_matrix(x))
    stop("'nrow' or 'ncol' cannot be specified when 'x' is a matrix")
  else {
    expr <- as.Constant(x)
    n <- length(expr)
    if(!missing(nrow))
      n <- nrow
    if(missing(ncol))
      ncol <- n
    expr*(base::diag(n)[1:n, 1:ncol])
  }
})
setMethod("diff", "Expression", function(x, lag = 1, differences = 1, ...) { Diff(x = x, lag = lag, k = differences, ...) })

setMethod("kronecker", signature(X = "Expression", Y = "ANY"), function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  if(FUN != "*" || make.dimnames)
    stop("Unimplemented")
  Kron(X, Y)
})
setMethod("kronecker", signature(X = "ANY", Y = "Expression"), function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  if(FUN != "*" || make.dimnames)
    stop("Unimplemented")
  Kron(X, Y)
})
setMethod("%x%", signature(X = "Expression", Y = "ANY"), function(X, Y) { Kron(lh_exp = X, rh_exp = Y) })
setMethod("%x%", signature(X = "ANY", Y = "Expression"), function(X, Y) { Kron(lh_exp = X, rh_exp = Y) })

setMethod("cbind2", signature(x = "Expression", y = "ANY"), function(x, y, ...) { HStack(x, y) })
setMethod("cbind2", signature(x = "ANY", y = "Expression"), function(x, y, ...) { HStack(x, y) })
setMethod("rbind2", signature(x = "Expression", y = "ANY"), function(x, y, ...) { VStack(x, y) })
setMethod("rbind2", signature(x = "ANY", y = "Expression"), function(x, y, ...) { VStack(x, y) })
