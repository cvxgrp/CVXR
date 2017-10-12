# =========================
#     Scalar functions
# =========================
geo_mean <- GeoMean

#'
#' Harmonic Mean
#'
#' The harmonic mean \eqn{\left(\frac{1}{n} \sum_{i=1}^n x_i^{-1}\right)^{-1}}. For a matrix, the function is applied over all entries.
#'
#' @param x An \S4class{Expression}, vector, or matrix.
#' @return An \S4class{Expression} representing the harmonic mean of the input.
#' @aliases harmonic_mean
#' @export
harmonic_mean <- HarmonicMean

#'
#' Maximum Eigenvalue
#' 
#' The maximum eigenvalue of a matrix \eqn{\lambda_{\max}(A)}.
#' 
#' @param A An \S4class{Expression} or matrix.
#' @return An \S4class{Expression} representing the maximum eigenvalue of the input.
#' @aliases lambda_max
#' @export
lambda_max <- LambdaMax

#'
#' Minimum Eigenvalue
#' 
#' The minimum eigenvalue of a matrix \eqn{\lambda_{\min}(A)}.
#' 
#' @param A An \S4class{Expression} or matrix.
#' @return An \S4class{Expression} representing the minimum eigenvalue of the input.
#' @aliases lambda_min
#' @export
lambda_min <- LambdaMin

#'
#' Sum of Largest Eigenvalues
#' 
#' The sum of the largest \eqn{k} eigenvalues of a matrix.
#' 
#' @param A An \S4class{Expression} or matrix.
#' @param k The number of eigenvalues to sum over.
#' @return An \S4class{Expression} representing the sum of the largest \code{k} eigenvalues of the input.
#' @aliases lambda_sum_largest
#' @export
lambda_sum_largest <- LambdaSumLargest

#'
#' Sum of Smallest Eigenvalues
#' 
#' The sum of the smallest \eqn{k} eigenvalues of a matrix.
#' 
#' @param A An \S4class{Expression} or matrix.
#' @param k The number of eigenvalues to sum over.
#' @return An \S4class{Expression} representing the sum of the smallest \code{k} eigenvalues of the input.
#' @aliases lambda_sum_smallest
#' @export
lambda_sum_smallest <- LambdaSumSmallest

#'
#' Log-Determinant
#' 
#' The natural logarithm of the determinant of a matrix \eqn{\log\det(A)}.
#' 
#' @param A An \S4class{Expression} or matrix.
#' @return An \S4class{Expression} representing the log-determinant of the input.
#' @aliases log_det
#' @export
log_det <- LogDet

#'
#' Log-Sum-Exponential
#' 
#' The natural logarithm of the sum of the elementwise exponential \eqn{\log\sum_{i=1}^n e^{x_i}}.
#' 
#' @param x An \S4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \S4class{Expression} representing the log-sum-exponential of the input.
#' @aliases log_sum_exp
#' @export
log_sum_exp <- LogSumExp

#'
#' Matrix Fraction
#' 
#' The matrix fraction function \eqn{tr(X^T P^{-1} X)}.
#' 
#' @param X An \S4class{Expression} or matrix. Must have the same number of rows as \code{P}.
#' @param P An \S4class{Expression} or matrix. Must be an invertible square matrix.
#' @return An \S4class{Expression} representing the matrix fraction evaluated at the input.
#' @aliases matrix_frac
#' @export
matrix_frac <- MatrixFrac

#'
#' Maximum
#' 
#' The maximum of an expression.
#' 
#' @param x An \S4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \S4class{Expression} representing the maximum of the input.
#' @aliases max_entries
#' @export
max_entries <- MaxEntries

#'
#' Minimum
#' 
#' The minimum of an expression.
#' 
#' @param x An \S4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \S4class{Expression} representing the minimum of the input.
#' @aliases min_entries
#' @export
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
#' @return An \S4class{Expression} representing the \eqn{l_{p,q}} norm of the input.
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
#' @return An \S4class{Expression} representing the 1-norm of the input.
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
#' @return An \S4class{Expression} representing the Euclidean norm of the input.
#' @aliases norm2
norm2 <- Norm2

#'
#' Infinity-Norm
#' 
#' The \eqn{\infty}}-norm \eqn{\|x\|_{\infty} = \max_{i=1,\ldots,n} |x_i|}.
#' 
#' @param x An \S4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \S4class{Expression} representing the infinity-norm of the input.
#' @aliases norm_inf
#' @export
norm_inf <- NormInf

#'
#' Nuclear Norm
#' 
#' The nuclear norm, i.e. sum of the singular values of a matrix.
#' 
#' @param A An \S4class{Expression} or matrix.
#' @return An \S4class{Expression} representing the nuclear norm of the input.
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
#' @return An \S4class{Expression} representing the quadratic form evaluated at the input.
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
#' @return An \S4class{Expression} representing the quadratic over linear function value evaluated at the input.
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
#' @return An \S4class{Expression} representing the sum of the entries of the input.
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
#' @return An \S4class{Expression} representing the sum of the largest k values of the input.
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
#' @return An \S4class{Expression} representing the sum of the smallest k values of the input.
#' @aliases sum_smallest
#' @export
sum_smallest <- SumSmallest

#'
#' Sum of Squares
#' 
#' The sum of the squared entries in a vector or matrix.
#' 
#' @param expr An \S4class{Expression}, vector, or matrix.
#' @return An \S4class{Expression} representing the sum of squares of the input.
#' @aliases sum_squares
#' @export
sum_squares <- SumSquares

#'
#' Matrix Trace
#'
#' The sum of the diagonal entries in a matrix.
#'
#' @param expr An \S4class{Expression} or matrix.
#' @return An \S4class{Expression} representing the trace of the input.
#' @aliases matrix_trace
matrix_trace <- Trace

#'
#' Total Variation
#' 
#' The total variation of a vector, matrix, or list of matrices. Uses L1 norm of discrete gradients for vectors and L2 norm of discrete gradients for matrices.
#' 
#' @param value An \S4class{Expression}, vector, or matrix.
#' @param ... \S4class{Expression} objects or numeric constants that extend the third dimension of value.
#' @return An \S4class{Expression} representing the total variation of the input.
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

# =========================
# Elementwise functions
# =========================

#'
#' Entropy Function
#'
#' The elementwise entropy function \eqn{-xlog(x)}.
#'
#' @param x An \S4class{Expression}, vector, or matrix.
#' @return An \S4class{Expression} representing the entropy of the input.
#' @aliases entr
#' @export
entr <- Entr
huber <- Huber

#'
#' Reciprocal Function
#'
#' The elementwise reciprocal function \eqn{\frac{1}{x}}
#'
#' @param x An \S4class{Expression}, vector, or matrix.
#' @return An \S4class{Expression} representing the reciprocal of the input.
#' @aliases inv_pos
#' @export 
inv_pos <- InvPos

#'
#' Kullback-Leibler Divergence
#'
#' The elementwise Kullback-Leibler divergence \eqn{x\log(x/y) - x + y}.
#'
#' @param x An \S4class{Expression}, vector, or matrix.
#' @param y An \S4class{Expression}, vector, or matrix.
#' @return An \S4class{Expression} representing the KL-divergence of the input.
#' @aliases kl_div
#' @export
kl_div <- KLDiv

#'
#' Logistic Function
#'
#' The elementwise logistic function \eqn{\log(1 + e^x)}.
#' This is a special case of log(sum(exp)) that evaluates to a vector rather than to a scalar, which is useful for logistic regression.
#'
#' @param x An \S4class{Expression}, vector, or matrix.
#' @return An \S4class{Expression} representing the logistic function evaluated at the input.
#' @aliases logistic
#' @export
logistic <- Logistic

#'
#' Elementwise Maximum
#'
#' The elementwise maximum.
#'
#' @param arg1 An \S4class{Expression}, vector, or matrix.
#' @param arg2 An \S4class{Expression}, vector, or matrix.
#' @param ... Additional \S4class{Expression} objects, vectors, or matrices.
#' @return An \S4class{Expression} 
#' @aliases max_elemwise
#' @export
max_elemwise <- MaxElemwise

#'
#' Elementwise Minimum
#'
#' The elementwise minimum.
#'
#' @param arg1 An \S4class{Expression}, vector, or matrix.
#' @param arg2 An \S4class{Expression}, vector, or matrix.
#' @param ... Additional \S4class{Expression} objects, vectors, or matrices.
#' @aliases min_elemwise
#' @export
min_elemwise <- MinElemwise

#'
#' Elementwise Multiplication
#'
#' The elementwise multiplication of two expressions. The first expression must be constant.
#'
#' @param lh_const A constant \S4class{Expression}, vector, or matrix representing the left-hand value.
#' @param rh_exp An \S4class{Expression}, vector, or matrix representing the right-hand value.
#' @aliases mul_elemwise
#' @export
mul_elemwise <- MulElemwise

#'
#' Elementwise Negative
#'
#' The elementwise absolute negative portion of an expression \eqn{-\min(x_i,0)}. This is equivalent to \code{-min_elemwise(x,0)}.
#'
#' @param x An \S4class{Expression}, vector, or matrix.
#' @return An \S4class{Expression} representing the negative portion of the input.
#' @aliases neg
#' @export
neg <- Neg

#'
#' Elementwise Positive
#'
#' The elementwise positive portion of an expression \eqn{\max(x_i,0)}. This is equivalent to \code{max_elemwise(x,0)}.
#'
#' @param x An \S4class{Expression}, vector, or matrix.
#' @return An \S4class{Expression} representing the positive portion of the input.
#' @aliases pos
#' @export
pos <- Pos
power <- Power

#'
#' Scalene Function
#'
#' The elementwise weighted sum of the positive and negative portions of an expression \eqn{\alpha\max(x_i,0) - \beta\min(x_i,0)}.
#' This is equivalent to \code{alpha*pos(x) + beta*neg(x)}.
#' 
#' @param x An \S4class{Expression}, vector, or matrix.
#' @param alpha The weight on the positive portion of \code{x}.
#' @param beta The weight on othe negative portion of \code{x}.
#' @return An \S4class{Expression} representing the scalene function evaluated at the input.
#' @aliases scalene
#' @export
scalene <- Scalene

#'
#' Square Function
#'
#' The elementwise square function \eqn{x^2}. This is equivalent to \code{power(x,2)}.
#' 
#' @param x An \S4class{Expression}, vector, or matrix.
#' @return An \S4class{Expression} representing the square of the input.
#' @aliases square
#' @export
square <- Square

setMethod("abs", "Expression", function(x) { Abs(x = x) })
setMethod("exp", "Expression", function(x) { Exp(x = x) })
setMethod("log", "Expression", function(x, base = exp(1)) { Log(x = x)/log(base) })
setMethod("log10", "Expression", function(x) { log(x, base = 10) })
setMethod("log2", "Expression", function(x) { log(x, base = 2) })
setMethod("log1p", "Expression", function(x) { Log1p(x = x) })
setMethod("sqrt", "Expression", function(x) { Sqrt(x = x) })

# =========================
# Matrix/vector operations
# =========================
bmat <- Bmat

#'
#' Discrete Convolution
#'
#' The 1-D discrete convolution of two vectors.
#'
#' @param lh_exp An \S4class{Expression} or vector representing the left-hand value.
#' @param rh_exp An \S4class{Expression} or vector representing the right-hand value.
#' @return An \S4class{Expression} representing the convolution of the input.
#' @aliases conv
#' @export
conv <- Conv
cum_sum <- CumSum

#'
#' Horizontal Concatenation
#'
#' The horizontal concatenation of expressions. This is equivalent to \code{cbind} when applied to objects with the same number of rows.
#' 
#' @param ... \S4class{Expression} objects, vectors, or matrices. All arguments must have the same number of rows.
#' @return An \S4class{Expression} representing the concatenated inputs.
#' @aliases hstack
#' @export
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

#'
#' Vectorization
#'
#' Flattens a matrix into a vector in column-major order.
#' 
#' @param X An \S4class{Expression} or matrix.
#' @return An \S4class{Expression} representing the vectorized matrix.
#' @aliases vec
#' @export
vec <- Vec

#'
#' Vertical Concatenation
#'
#' The vertical concatenation of expressions. This is equivalent to \code{rbind} when applied to objects with the same number of columns.
#' 
#' @param ... \S4class{Expression} objects, vectors, or matrices. All arguments must have the same number of columns.
#' @return An \S4class{Expression} representing the concatenated inputs.
#' @aliases vstack
#' @export
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
