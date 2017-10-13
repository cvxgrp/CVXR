# =========================
#     Scalar functions
# =========================
#'
#' Geometric Mean
#'
#' The (weighted) geometric mean of vector \eqn{x} with optional powers given by \eqn{p}.
#'
#' \deqn{\left(x_1^{p_1} \cdots x_n^{p_n} \right)^{\frac{1}{\mathbf{1}^Tp}}}
#'
#' The geometric mean includes an implicit constraint that \eqn{x_i \geq 0} whenever \eqn{p_i > 0}. If \eqn{p_i = 0, x_i} will be unconstrained.
#' The only exception to this rule occurs when \eqn{p} has exactly one nonzero element, say \eqn{p_i}, in which case \code{geo_mean(x,p)} is equivalent to \eqn{x_i} (without the nonnegativity constraint).
#' A specific case of this is when \eqn{x \in \mathbf{R}^1}.
#'
#' @param x An \linkS4class{Expression} or vector.
#' @param p (Optional) A vector of weights for the weighted geometric mean. Defaults to a vector of ones, giving the \strong{unweighted} geometric mean \eqn{x_1^{1/n} \cdots x_n^{1/n}}.
#' @return An \linkS4class{Expression} representing the geometric mean of the input.
#' @examples 
#' x <- Variable(2)
#' cost <- geo_mean(x)
#' prob <- Problem(Maximize(cost), list(sum(x) <= 1))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' 
#' prob <- Problem(Maximize(geo_mean(x, p)), list(p_norm(x) <= 1))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' @docType methods
#' @rdname geo_mean
#' @export
geo_mean <- GeoMean

#'
#' Harmonic Mean
#'
#' The harmonic mean, \eqn{\left(\frac{1}{n} \sum_{i=1}^n x_i^{-1}\right)^{-1}}. For a matrix, the function is applied over all entries.
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @return An \linkS4class{Expression} representing the harmonic mean of the input.
#' @examples 
#' x <- Variable(3)
#' val <- c(1,2,3)
#' prob <- Problem(Minimize(harmonic_mean(x)), list(x == val))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @rdname harmonic_mean
#' @export
harmonic_mean <- HarmonicMean

#'
#' Maximum Eigenvalue
#'
#' The maximum eigenvalue of a matrix, \eqn{\lambda_{\max}(A)}.
#' @param A An \linkS4class{Expression} or matrix.
#' @return An \linkS4class{Expression} representing the maximum eigenvalue of the input.
#' @examples 
#' A <- Variable(2, 2, name = "A")
#' p <- Problem(Minimize(lambda_max(A)), list(A >= 2))
#' result <- solve(p)
#' result$value
#' result$getValue(A)
#' 
#' obj <- Maximize(A[2,1] - A[1,2])
#' p <- Problem(obj, list(lambda_max(A) <= 100, A[1,1] == 2, A[2,2] == 2, A[2,1] == 2))
#' result <- solve(p)
#' result$value
#' result$getValue(A)
#' @docType methods
#' @rdname lambda_max
#' @export
lambda_max <- LambdaMax

#'
#' Minimum Eigenvalue
#'
#' The minimum eigenvalue of a matrix, \eqn{\lambda_{\min}(A)}.
#' 
#' @param A An \linkS4class{Expression} or matrix.
#' @return An \linkS4class{Expression} representing the minimum eigenvalue of the input.
#' @examples 
#' A <- Variable(2, 2, name = "A")
#' val <- cbind(c(5,7), c(7,-3))
#' p <- Problem(Minimize(lambda_min(A)), list(A == val))
#' result <- solve(p)
#' result$value
#' result$getValue(A)
#' @docType methods
#' @rdname lambda_min
#' @export
lambda_min <- LambdaMin

#'
#' Sum of Largest Eigenvalues
#'
#' The sum of the largest \eqn{k} eigenvalues of a matrix.
#' 
#' @param A An \linkS4class{Expression} or matrix.
#' @param k The number of eigenvalues to sum over.
#' @return An \linkS4class{Expression} representing the sum of the largest \code{k} eigenvalues of the input.
#' @docType methods
#' @rdname lambda_sum_largest
#' @export
lambda_sum_largest <- LambdaSumLargest

#'
#' Sum of Smallest Eigenvalues
#'
#' The sum of the smallest \eqn{k} eigenvalues of a matrix.
#' 
#' @param A An \linkS4class{Expression} or matrix.
#' @param k The number of eigenvalues to sum over.
#' @return An \linkS4class{Expression} representing the sum of the smallest \code{k} eigenvalues of the input.
#' @docType methods
#' @rdname lambda_sum_smallest
#' @export
lambda_sum_smallest <- LambdaSumSmallest

#'
#' Log-Determinant
#'
#' The natural logarithm of the determinant of a matrix, \eqn{\log\det(A)}.
#' 
#' @param A An \linkS4class{Expression} or matrix.
#' @return An \linkS4class{Expression} representing the log-determinant of the input.
#' @docType methods
#' @rdname log_det
#' @export
log_det <- LogDet

#'
#' Log-Sum-Exponential
#'
#' The natural logarithm of the sum of the elementwise exponential, \eqn{\log\sum_{i=1}^n e^{x_i}}.
#' 
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \linkS4class{Expression} representing the log-sum-exponential of the input.
#' @docType methods
#' @rdname log_sum_exp
#' @export
log_sum_exp <- LogSumExp

#'
#' Matrix Fraction
#'
#' \eqn{tr(X^T P^{-1} X)}.
#' 
#' @param X An \linkS4class{Expression} or matrix. Must have the same number of rows as \code{P}.
#' @param P An \linkS4class{Expression} or matrix. Must be an invertible square matrix.
#' @return An \linkS4class{Expression} representing the matrix fraction evaluated at the input.
#' @docType methods
#' @rdname matrix_frac
#' @export
matrix_frac <- MatrixFrac

#'
#' Maximum
#'
#' The maximum of an expression.
#' 
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \linkS4class{Expression} representing the maximum of the input.
#' @docType methods
#' @rdname max_entries
#' @export
max_entries <- MaxEntries

#'
#' Minimum
#'
#' The minimum of an expression.
#' 
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \linkS4class{Expression} representing the minimum of the input.
#' @docType methods
#' @rdname min_entries
#' @export
min_entries <- MinEntries

#'
#' Mixed Norm
#'
#' \eqn{l_{p,q}(x) = \left(\sum_{i=1}^n (\sum_{j=1}^m |x_{i,j}|)^{q/p}\right)^{1/q}}.
#'
#' @param X An \linkS4class{Expression}, vector, or matrix.
#' @param p The type of inner norm.
#' @param q The type of outer norm.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \linkS4class{Expression} representing the \eqn{l_{p,q}} norm of the input.
#' @docType methods
#' @rdname mixed_norm
#' @export
mixed_norm <- MixedNorm

#'
#' 1-Norm
#'
#' \eqn{\|x\|_1 = \sum_{i=1}^n |x_i|}.
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \linkS4class{Expression} representing the 1-norm of the input.
#' @examples
#' a <- Variable(name = "a")
#' p <- Problem(Minimize(norm1(a)), list(a <= -2))
#' result <- solve(p)
#' result$value
#' result$getValue(a)
#' 
#' p <- Problem(Maximize(-norm1(a)), list(a <= -2))
#' result <- solve(p)
#' result$value
#' result$getValue(a)
#' 
#' x <- Variable(2, name = "x")
#' z <- Variable(2, name = "z")
#' p <- Problem(Minimize(norm1(x - z) + 5), list(x >= c(2,3), z <= c(-1,-4)))
#' result <- solve(p)
#' result$value
#' result$getValue(x[1] - z[1])
#' @docType methods
#' @rdname norm1
#' @export
norm1 <- Norm1

#'
#' Euclidean Norm
#'
#' \eqn{\|x\|_2 = \left(\sum_{i=1}^n x_i^2\right)^{1/2}}.
#' 
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \linkS4class{Expression} representing the Euclidean norm of the input.
#' @examples
#' a <- Variable(name = "a")
#' p <- Problem(Minimize(norm2(a)), list(a <= -2))
#' result <- solve(p)
#' result$value
#' result$getValue(a)
#'
#' p <- Problem(Maximize(-norm2(a)), list(a <= -2))
#' result <- solve(p)
#' result$value
#' result$getValue(a)
#'
#' x <- Variable(2, name = "x")
#' z <- Variable(2, name = "z")
#' p <- Problem(Minimize(norm2(x - z) + 5), list(x >= c(2,3), z <= c(-1,-4)))
#' result <- solve(p)
#' result$value
#' result$getValue(x)
#' result$getValue(z)
#'
#' p <- Problem(Minimize(norm2(t(x - z)) + 5), list(x >= c(2,3), z <= c(-1,-4)))
#' result <- solve(p)
#' result$value
#' result$getValue(x)
#' result$getValue(z)
#' @docType methods
#' @rdname norm2
norm2 <- Norm2

#'
#' Infinity-Norm
#'
#' \eqn{\|x\|_{\infty} = \max_{i=1,\ldots,n} |x_i|}.
#' 
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \linkS4class{Expression} representing the infinity-norm of the input.
#' @examples
#' a <- Variable(name = "a")
#' b <- Variable(name = "b")
#' c <- Variable(name = "c")
#' 
#' p <- Problem(Minimize(norm_inf(a)), list(a >= 2))
#' result <- solve(p)
#' result$value
#' result$getValue(a)
#' 
#' p <- Problem(Minimize(3*norm_inf(a + 2*b) + c), list(a >= 2, b <= -1, c == 3))
#' result <- solve(p)
#' result$value
#' result$getValue(a + 2*b)
#' result$getValue(c)
#' 
#' p <- Problem(Maximize(-norm_inf(a)), list(a <= -2))
#' result <- solve(p)
#' result$value
#' result$getValue(a)
#'
#' x <- Variable(2, name = "x")
#' z <- Variable(2, name = "z")
#' p <- Problem(Minimize(norm_inf(x - z) + 5), list(x >= c(2,3), z <= c(-1,-4)))
#' result <- solve(p)
#' result$value
#' result$getValue(x[1] - z[1])
#' @docType methods
#' @rdname norm_inf
#' @export
norm_inf <- NormInf

#'
#' Nuclear Norm
#'
#' The nuclear norm, i.e. sum of the singular values of a matrix.
#' 
#' @param A An \linkS4class{Expression} or matrix.
#' @return An \linkS4class{Expression} representing the nuclear norm of the input.
#' @docType methods
#' @rdname norm_nuc
#' @export
norm_nuc <- NormNuc

#'
#' P-Norm
#'
#' The vector p-norm. If given a matrix variable, \code{p_norm} will treat it as a vector and compute the p-norm of the concatenated columns.
#'
#' For \eqn{p \geq 1}, the p-norm is given by \deqn{\|x\|_p = \left(\sum_{i=1}^n |x_i|^p\right)^{1/p}} with domain \eqn{x \in \mathbf{R}^n}.
#' For \eqn{p < 1, p \neq 0}, the p-norm is given by \deqn{\|x\|_p = \left(\sum_{i=1}^n x_i^p\right)^{1/p}} with domain \eqn{x \in \mathbf{R}^n_+}.
#'
#' \itemize{
#'    \item Note that the "p-norm" is actually a \strong{norm} only when \eqn{p \geq 1} or \eqn{p = +\infty}. For these cases, it is convex.
#'    \item The expression is undefined when \eqn{p = 0}.
#'    \item Otherwise, when \eqn{p < 1}, the expression is concave, but not a true norm.
#' }
#' 
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param p A number greater than or equal to 1, or equal to positive infinity.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code[NA}.
#' @return An \linkS4class{Expression} representing the p-norm of the input.
#' @examples 
#' x <- Variable(3, name = "x")
#' prob <- Problem(Minimize(p_norm(x, 2)))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' 
#' prob <- Problem(Minimize(p_norm(x, Inf)))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#'
#' a <- c(1.0, 2, 3)
#' prob <- Problem(Minimize(p_norm(x, 1.6)), list(t(x) %*% a >= 1))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' 
#' prob <- Problem(Minimize(sum(abs(x - a))), list(p_norm(x,-1) >= 0))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' @docType methods
#' @rdname p_norm
#' @export
p_norm <- Pnorm

#'
#' Quadratic Form
#'
#' The quadratic form, \eqn{x^TPx}.
#'
#' @param x An \linkS4class{Expression} or vector.
#' @param P An \linkS4class{Expression} or matrix.
#' @return An \linkS4class{Expression} representing the quadratic form evaluated at the input.
#' @examples 
#' x <- Variable(2, name = "x")
#' P <- rbind(c(4,0), c(0,9))
#' p <- Problem(Minimize(quad_form(x, P)), list(x >= 1))
#' result <- solve(p)
#' result$value
#' result$getValue(x)
#'
#' A <- Variable(2, 2, name = "A")
#' c <- c(1,2)
#' p <- Problem(Minimize(quad_form(c, A)), list(A >= 1))
#' result <- solve(p)
#' result$value
#' result$getValue(A)
#' @docType methods
#' @rdname quad_form
#' @export
quad_form <- QuadForm

#'
#' Quadratic over Linear
#'
#' \eqn{\sum_{i,j} X_{i,j}^2/y}.
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param y A scalar \linkS4class{Expression} or numeric constant.
#' @return An \linkS4class{Expression} representing the quadratic over linear function value evaluated at the input.
#' @docType methods
#' @rdname quad_over_lin
#' @export
quad_over_lin <- QuadOverLin

#'
#' Sum of Entries
#'
#' The sum of entries in a vector or matrix.
#'
#' @param expr An \linkS4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code[NA}.
#' @return An \linkS4class{Expression} representing the sum of the entries of the input.
#' @docType methods
#' @rdname sum_entries
#' @export
sum_entries <- SumEntries

#'
#' Sum of Largest Values
#'
#' The sum of the largest \eqn{k} values of a vector or matrix.
#' 
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param k The number of largest values to sum over.
#' @return An \linkS4class{Expression} representing the sum of the largest \code{k} values of the input.
#' @docType methods
#' @rdname sum_largest
#' @export
sum_largest <- SumLargest

#'
#' Sum of Smallest Values
#'
#' The sum of the smallest k values of a vector or matrix.
#' 
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param k The number of smallest values to sum over.
#' @return An \linkS4class{Expression} representing the sum of the smallest k values of the input.
#' @docType methods
#' @rdname sum_smallest
#' @export
sum_smallest <- SumSmallest

#'
#' Sum of Squares
#'
#' The sum of the squared entries in a vector or matrix.
#' 
#' @param expr An \linkS4class{Expression}, vector, or matrix.
#' @return An \linkS4class{Expression} representing the sum of squares of the input.
#' @docType methods
#' @rdname sum_squares
#' @export
sum_squares <- SumSquares

#'
#' Matrix Trace
#'
#' The sum of the diagonal entries in a matrix.
#'
#' @param expr An \linkS4class{Expression} or matrix.
#' @return An \linkS4class{Expression} representing the trace of the input.
#' @docType methods
#' @rdname matrix_trace
matrix_trace <- Trace

#'
#' Total Variation
#'
#' The total variation of a vector, matrix, or list of matrices. Uses L1 norm of discrete gradients for vectors and L2 norm of discrete gradients for matrices.
#' 
#' @param value An \linkS4class{Expression}, vector, or matrix.
#' @param ... (Optional) \linkS4class{Expression} objects or numeric constants that extend the third dimension of value.
#' @return An \linkS4class{Expression} representing the total variation of the input.
#' @docType methods
#' @rdname tv
#' @export
tv <- TotalVariation

#' @docType methods
#' @rdname max_entries
#' @export
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

#' @docType methods
#' @rdname min_entries
#' @export
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

#'
#' Matrix Norm
#'
#' The matrix norm, which can be the 1-norm ("1"), infinity-norm ("I"), Frobenius norm ("F"), maximum modulus of all the entries ("M"), or the spectral norm ("2"), as determined by the value of type.
#' 
#' @param x An \S4class{Expression}.
#' @type A character indicating the type of norm desired.
#' \itemize{
#'    \item "O", "o" or "1" specifies the 1-norm (maximum absolute column sum).
#'    \item "I" or "i" specifies the infinity-norm (maximum absolute row sum).
#'    \item "F" or "f" specifies the Frobenius norm (Euclidean norm of the vectorized \code{x}).
#'    \item "M" or "m" specifies the maximum modulus of all the elements in \code{x}.
#'    \item "2" specifies the spectral norm, which is the largest singular value of \code{x}.
#' }
#' @return An \S4class{Expression} representing the norm of the input.
#' @seealso The \code{\link{p_norm}} function calculates the vector p-norm.
#' @docType methods
#' @rdname norm
#' @export
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

#' @docType methods
#' @rdname sum_entries
#' @export
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

#'
#' Arithmetic Mean
#'
#' The arithmetic mean of an expression.
#' 
#' @param x An \linkS4class{Expression}.
#' @return An \linkS4class{Expression} representing the mean of the input.
#' @docType methods
#' @rdname mean
#' @export
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
#' The elementwise entropy function, \eqn{-xlog(x)}.
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @return An \linkS4class{Expression} representing the entropy of the input.
#' @docType methods
#' @rdname entr
#' @export
entr <- Entr

#'
#' Huber Function
#'
#' The elementwise Huber function,
#' \deqn{\mbox{Huber}(x, M) = \begin{cases}
#'       2M|x|-M^2 & \mbox{for } |x| \geq |M| \\
#'       |x|^2 & \mbox{for } |x| \leq M
#' \end{cases}}
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param M (Optional) A positive scalar value representing the threshold. Defaults to 1.
#' @return An \linkS4class{Expression} representing the Huber function evaluated at the input.
#' @docType methods
#' @rdname huber
#' @export
huber <- Huber

#'
#' Reciprocal Function
#'
#' The elementwise reciprocal function, \eqn{\frac{1}{x}}
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @return An \linkS4class{Expression} representing the reciprocal of the input.
#' @docType methods
#' @rdname inv_pos
#' @export
inv_pos <- InvPos

#'
#' Kullback-Leibler Divergence
#'
#' The elementwise Kullback-Leibler divergence, \eqn{x\log(x/y) - x + y}.
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param y An \linkS4class{Expression}, vector, or matrix.
#' @return An \linkS4class{Expression} representing the KL-divergence of the input.
#' @docType methods
#' @rdname kl_div
#' @export
kl_div <- KLDiv

#'
#' Logistic Function
#'
#' The elementwise logistic function, \eqn{\log(1 + e^x)}.
#' This is a special case of log(sum(exp)) that evaluates to a vector rather than to a scalar, which is useful for logistic regression.
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @return An \linkS4class{Expression} representing the logistic function evaluated at the input.
#' @docType methods
#' @rdname logistic
#' @export
logistic <- Logistic

#'
#' Elementwise Maximum
#'
#' The elementwise maximum.
#'
#' @param arg1 An \linkS4class{Expression}, vector, or matrix.
#' @param arg2 An \linkS4class{Expression}, vector, or matrix.
#' @param ... Additional \linkS4class{Expression} objects, vectors, or matrices.
#' @return An \linkS4class{Expression} representing the elementwise maximum of the inputs.
#' @docType methods
#' @rdname max_elemwise
#' @export
max_elemwise <- MaxElemwise

#'
#' Elementwise Minimum
#'
#' The elementwise minimum.
#'
#' @param arg1 An \linkS4class{Expression}, vector, or matrix.
#' @param arg2 An \linkS4class{Expression}, vector, or matrix.
#' @param ... Additional \linkS4class{Expression} objects, vectors, or matrices.
#' @return An \linkS4class{Expression} representing the elementwise minimum of the inputs.
#' @docType methods
#' @rdname min_elemwise
#' @export
min_elemwise <- MinElemwise

#'
#' Elementwise Multiplication
#'
#' The elementwise product of two expressions. The first expression must be constant.
#'
#' @param lh_const A constant \linkS4class{Expression}, vector, or matrix representing the left-hand value.
#' @param rh_exp An \linkS4class{Expression}, vector, or matrix representing the right-hand value.
#' @return An \linkS4class{Expression} representing the elementwise product of the inputs.
#' @examples 
#' A <- Variable(2, 2, name = "A")
#' c <- cbind(c(1,-1), c(2,-2))
#' expr <- mul_elemwise(c, A)
#' obj <- Minimize(norm_inf(expr))
#' p <- Problem(obj, list(A == 5))
#' result <- solve(p)
#' result$value
#' result$getValue(expr)
#' @docType methods
#' @rdname mul_elemwise
#' @export
mul_elemwise <- MulElemwise

#'
#' Elementwise Negative
#'
#' The elementwise absolute negative portion of an expression, \eqn{-\min(x_i,0)}. This is equivalent to \code{-min_elemwise(x,0)}.
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @return An \linkS4class{Expression} representing the negative portion of the input.
#' @docType methods
#' @rdname neg
#' @export
neg <- Neg

#'
#' Elementwise Positive
#'
#' The elementwise positive portion of an expression, \eqn{\max(x_i,0)}. This is equivalent to \code{max_elemwise(x,0)}.
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @return An \linkS4class{Expression} representing the positive portion of the input.
#' @docType methods
#' @rdname pos
#' @export
pos <- Pos

#'
#' Elementwise Power
#'
#' Raises each element of the input to the power \eqn{p}. 
#' If \code{expr} is a CVXR expression, then \code{expr^p} is equivalent to \code{power(expr,p)}.
#'
#' For \eqn{p = 0} and \eqn{f(x) = 1}, this function is constant and positive.
#' For \eqn{p = 1} and \eqn{f(x) = x}, this function is affine, increasing, and the same sign as \eqn{x}.
#' For \eqn{p = 2,4,8,\ldots} and \eqn{f(x) = |x|^p}, this function is convex, positive, with signed monotonicity.
#' For \eqn{p < 0} and \deqn{f(x) = \begin{cases} x^p & x > 0 \\ +\infty & x \leq 0 \end{cases}}, this function is convex, decreasing, and positive.
#' For \eqn{0 < p < 1} and \deqn{f(x) = \begin{cases} x^p & x \geq 0 \\ -\infty & x < 0 \end{cases}}, this function is concave, increasing, and positive.
#' For \eqn{p > 1, p \neq 2,4,8,\ldots} and \deqn{f(x) = \begin{cases} x^p & x \geq 0 \\ +\infty & x < 0 \end{cases}}, this function is convex, increasing, and positive.
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param p A scalar value indicating the exponential power.
#' @param max_denom The maximum denominator considered in forming a rational approximation of \code{p}.
#' @examples 
#' x <- Variable()
#' prob <- Problem(Minimize(power(x, 1.7) + power(x, -2.3) - power(x, 0.45)))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' @docType methods
#' @rdname power
#' @export
power <- Power

#'
#' Scalene Function
#'
#' The elementwise weighted sum of the positive and negative portions of an expression, \eqn{\alpha\max(x_i,0) - \beta\min(x_i,0)}.
#' This is equivalent to \code{alpha*pos(x) + beta*neg(x)}.
#' 
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param alpha The weight on the positive portion of \code{x}.
#' @param beta The weight on othe negative portion of \code{x}.
#' @return An \linkS4class{Expression} representing the scalene function evaluated at the input.
#' @docType methods
#' @rdname scalene
#' @export
scalene <- Scalene

#'
#' Square Function
#'
#' The elementwise square function. This is equivalent to \code{power(x,2)}.
#' 
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @return An \linkS4class{Expression} representing the square of the input.
#' @docType methods
#' @rdname square
#' @export
square <- Square

#'
#' Absolute Value
#'
#' The elementwise absolute value.
#'
#' @param x An \linkS4class{Expression}.
#' @return An \linkS4class{Expression} representing the absolute value of the input.
#' @examples
#' A <- Variable(2, 2, name = "A")
#' p <- Problem(Minimize(sum(abs(A))), list(-2 >= A))
#' result <- solve(p)
#' result$value
#' result$getValue(A)
#' @docType methods
#' @rdname abs
#' @export
setMethod("abs", "Expression", function(x) { Abs(x = x) })

#'
#' Natural Exponential
#'
#' The elementwise natural exponential.
#'
#' @param x An \linkS4class{Expression}.
#' @return An \linkS4class{Expression} representing the natural exponential of the input.
#' @docType methods
#' @rdname exp
#' @export
setMethod("exp", "Expression", function(x) { Exp(x = x) })

#'
#' Logarithms
#'
#' The elementwise logarithm.
#' \code{log} computes the logarithm, by default the natural logarithm, \code{log10} computes the common (i.e., base 10) logarithm, and \code{log2} computes the binary (i.e., base 2) logarithms. The general form \code{log(x, base)} computes logarithms with base \code{base}.
#' \code{log1p} computes elementwise the function \eqn{\log(1+x)}.
#'
#' @param x An \linkS4class{Expression}.
#' @param base (Optional) A positive number that is the base with respect to which the logarithm is computed. Defaults to \eqn{e}.
#' @return An \linkS4class{Expression} representing the exponentiated input.
#' @docType methods
#' @rdname log
#' @export
setMethod("log", "Expression", function(x, base = exp(1)) { Log(x = x)/log(base) })

#' @docType methods
#' @rdname log
#' @export
setMethod("log10", "Expression", function(x) { log(x, base = 10) })

#' @docType methods
#' @rdname log
#' @export
setMethod("log2", "Expression", function(x) { log(x, base = 2) })

#' @docType methods
#' @rdname log
#' @export
setMethod("log1p", "Expression", function(x) { Log1p(x = x) })
log1p <- Log1p

#'
#' Square Root
#'
#' The elementwise square root.
#' 
#' @param x An \linkS4class{Expression}.
#' @return An \linkS4class{Expression} representing the square root of the input.
#' @rdname sqrt
#' @export
setMethod("sqrt", "Expression", function(x) { Sqrt(x = x) })

# =========================
# Matrix/vector operations
# =========================
#'
#' Block Matrix
#'
#' Constructs a block matrix from a list of lists. Each internal list is stacked horizontally, and the internal lists are stacked vertically.
#' 
#' @param block_lists A list of lists containing \linkS4class{Expression} objects, matrices, or vectors, which represent the blocks of the block matrix.
#' @return An \linkS4class{Expression} representing the block matrix.
#' @docType methods
#' @rdname bmat
#' @export
bmat <- Bmat

#'
#' Discrete Convolution
#'
#' The 1-D discrete convolution of two vectors.
#'
#' @param lh_exp An \linkS4class{Expression} or vector representing the left-hand value.
#' @param rh_exp An \linkS4class{Expression} or vector representing the right-hand value.
#' @return An \linkS4class{Expression} representing the convolution of the input.
#' @examples 
#' x <- Variable(5)
#' h <- matrix(rnorm(2), nrow = 2, ncol = 1)
#' prob <- Problem(Minimize(sum(conv(h, x))))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' @docType methods
#' @rdname conv
#' @export
conv <- Conv

#'
#' Horizontal Concatenation
#'
#' The horizontal concatenation of expressions.
#' This is equivalent to \code{cbind} when applied to objects with the same number of rows.
#' 
#' @param ... \linkS4class{Expression} objects, vectors, or matrices. All arguments must have the same number of rows.
#' @return An \linkS4class{Expression} representing the concatenated inputs.
#' @examples 
#' x <- Variable(2, name = "x")
#' y <- Variable(3, name = "y")
#' c <- matrix(1, nrow = 1, ncol = 5)
#' p <- Problem(Minimize(c %*% t(hstack(t(x), t(y)))), list(x == c(1,2), y == c(3,4,5)))
#' result <- solve(p)
#' result$value
#'
#' c <- matrix(1, nrow = 1, ncol = 4)
#' p <- Problem(Minimize(c %*% t(hstack(t(x), t(x)))), list(x == c(1,2)))
#' result <- solve(p)
#' result$value
#'
#' A <- Variable(2, 2, name = "A")
#' C <- Variable(3, 2, name = "C")
#' c <- matrix(1, nrow = 2, ncol = 2)
#' p <- Problem(Minimize(sum_entries(hstack(t(A), t(C)))), list(A >= 2*c, C == -2))
#' result <- solve(p)
#' result$value
#' result$getValue(A)
#'
#' D <- Variable(3,3)
#' expr <- hstack(C, D)
#' p <- Problem(Minimize(expr[1,2] + sum(hstack(expr, expr))), list(C >= 0, D >= 0, D[1,1] == 2, C[1,2] == 3))
#' result <- solve(p)
#' result$value
#' result$getValue(C)
#' result$getValue(D)
#' @docType methods
#' @rdname hstack
#' @export
hstack <- HStack

#'
#' Reshape an Expression
#'
#' This function vectorizes an expression, then unvectorizes it into a new shape. Entries are stored in column-major order.
#'
#' @param expr An \linkS4class{Expression}, vector, or matrix.
#' @param rows The new number of rows.
#' @param cols The new number of columns.
#' @return An \linkS4class{Expression} representing the reshaped input.
#' @examples 
#' x <- Variable(4)
#' mat <- cbind(c(1,-1), c(2,-2))
#' vec <- matrix(1:4)
#' expr <- reshape_expr(x,2,2)
#' obj <- Minimize(sum(mat %*% expr))
#' prob <- Problem(obj, list(x == vec))
#' result <- solve(prob)
#' result$value
#' 
#' A <- Variable(2, 2, name = "A")
#' c <- 1:4
#' expr <- reshape_expr(A,4,1)
#' obj <- Minimize(t(expr) %*% c)
#' constraints <- list(A == cbind(c(-1,-2), c(3,4)))
#' prob <- Problem(obj, constraints)
#' result <- solve(prob)
#' result$value
#' result$getValue(expr)
#' result$getValue(reshape_expr(expr,2,2))
#'
#' C <- Variable(3, 2, name = "C")
#' expr <- reshape_expr(C,2,3)
#' mat <- rbind(c(1,-1), c(2,-2))
#' C_mat <- rbind(c(1,4), c(2,5), c(3,6))
#' obj <- Minimize(sum(mat %*% expr))
#' prob <- Problem(obj, list(C == C_mat))
#' result <- solve(prob)
#' result$value
#' result$getValue(expr)
#' 
#' a <- Variable(name = "a")
#' c <- cbind(c(1,-1), c(2,-2))
#' expr <- reshape_expr(c * a,1,4)
#' obj <- Minimize(expr %*% (1:4))
#' prob <- Problem(obj, list(a == 2))
#' result <- solve(prob)
#' result$value
#' result$getValue(expr)
#' 
#' expr <- reshape_expr(c * a,4,1)
#' obj <- Minimize(t(expr) %*% (1:4))
#' prob <- Problem(obj, list(a == 2))
#' result <- solve(prob)
#' result$value
#' result$getValue(expr)
#' @docType methods
#' @rdname reshape_expr
#' @export
reshape_expr <- Reshape

#'
#' Maximum Singular Value
#'
#' The maximum singular value of a matrix.
#'
#' @param A An \linkS4class{Expression} or matrix.
#' @return An \linkS4class{Expression} representing the maximum singular value.
#' @docType methods
#' @rdname sigma_max
#' @export
sigma_max <- SigmaMax

#'
#' Upper Triangle of a Matrix
#'
#' The vectorized strictly upper triangular entries of a matrix.
#'
#' @param expr An \linkS4class{Expression} or matrix.
#' @return An \linkS4class{Expression} representing the upper triangle of the input.
#' @docType methods
#' @rdname upper_tri
#' @export
upper_tri <- UpperTri

#'
#' Vectorization of a Matrix
#'
#' Flattens a matrix into a vector in column-major order.
#' 
#' @param X An \linkS4class{Expression} or matrix.
#' @return An \linkS4class{Expression} representing the vectorized matrix.
#' @examples
#' A <- Variable(2, 2, name = "A")
#' c <- 1:4
#' expr <- vec(A)
#' obj <- Minimize(t(expr) %*% c)
#' constraints <- list(A == cbind(c(-1,-2), c(3,4)))
#' prob <- Problem(obj, constraints)
#' result <- solve(prob)
#' result$value
#' result$getValue(expr)
#' @docType methods
#' @rdname vec
#' @export
vec <- Vec

#'
#' Vertical Concatenation
#'
#' The vertical concatenation of expressions. This is equivalent to \code{rbind} when applied to objects with the same number of columns.
#' 
#' @param ... \linkS4class{Expression} objects, vectors, or matrices. All arguments must have the same number of columns.
#' @return An \linkS4class{Expression} representing the concatenated inputs.
#' @examples 
#' x <- Variable(2, name = "x")
#' y <- Variable(3, name = "y")
#' c <- matrix(1, nrow = 1, ncol = 5)
#' p <- Problem(Minimize(c %*% vstack(x, y)), list(x == c(1,2), y == c(3,4,5)))
#' result <- solve(p)
#' result$value
#'
#' c <- matrix(1, nrow = 1, ncol = 4)
#' p <- Problem(Minimize(c %*% vstack(x, x)), list(x == c(1,2)))
#' result <- solve(p)
#' result$value
#'
#' A <- Variable(2, 2, name = "A")
#' C <- Variable(3, 2, name = "C")
#' c <- matrix(1, nrow = 2, ncol = 2)
#' p <- Problem(Minimize(sum(vstack(A, C))), list(A >= 2*c, C == -2))
#' result <- solve(p)
#' result$value
#'
#' B <- Variable(2, 2, name = "B")
#' c <- matrix(1, nrow = 1, ncol = 2)
#' p <- Problem(Minimize(sum(vstack(c %*% A, c %*% B))), list(A >= 2, B == -2))
#' result <- solve(p)
#' result$value
#' @docType methods
#' @rdname vstack
#' @export
vstack <- VStack

#'
#' Cumulative Sum
#'
#' The cumulative sum, \eqn{\sum_{i=1}^k x_i}. Matrices are flattened into column-major order before the sum is taken.
#'
#' @param expr An \linkS4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @docType methods
#' @rdname cumsum
#' @export
setMethod("cumsum", "Expression", function(x) { CumSum(expr = Vec(x)) })   # Flatten matrix in column-major order to match R's behavior
cumsum <- CumSum

#'
#' Matrix Diagonal
#'
#' Extracts the diagonal from a matrix or makes a vector into a diagonal matrix.
#'
#' @param expr An \linkS4class{Expression}, vector, or square matrix.
#' @return An \linkS4class{Expression} representing the diagonal vector or matrix.
#' @examples 
#' C <- Variable(3,3)
#' obj <- Maximize(C[1,3])
#' constraints <- list(diag(C) == 1, C[1,2] == 0.6, C[2,3] == -0.3, C == Semidef(3))
#' prob <- Problem(obj, constraints)
#' result <- solve(prob)
#' result$value
#' result$getValue(C)
#' @docType methods
#' @rdbane diag
#' @export
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

#'
#' Lagged and Iterated Differences
#'
#' The lagged and iterated differences of a vector.
#' If \code{x} is length \code{n}, this function returns a length \eqn{n-k} vector of the \eqn{k}th order difference between the lagged terms.
#' \code{diff(x)} returns the vector of differences between adjacent elements in the vector, i.e. [x[2] - x[1], x[3] - x[2], ...].
#' \code{diff(x,1,2)} is the second-order differences vector, equivalently diff(diff(x)). \code{diff(x,1,0)} returns the vector x unchanged.
#' \code{diff(x,2)} returns the vector of differences [x[3] - x[1], x[4] - x[2], ...], equivalent to \code{x[(1+lag):n] - x[1:(n-lag)]}.
#' 
#' @param x An \linkS4class{Expression}.
#' @param lag An integer indicating which lag to use.
#' @param differences An integer indicating the order of the difference.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \linkS4class{Expression} representing the \code{k}th order difference.
#' @docType methods
#' @rdname diff
#' @export
setMethod("diff", "Expression", function(x, lag = 1, differences = 1, ...) { Diff(x = x, lag = lag, k = differences, ...) })

#'
#' Kronecker Product
#' 
#' The generalized kronecker product of two matrices.
#' 
#' @param X An \linkS4class{Expression} or matrix.
#' @param Y An \linkS4class{Expression} or matrix.
#' @return An \linkS4class{Expression} that represents the kronecker product.
#' @docType methods
#' @rdname kronecker
#' @export
setMethod("kronecker", signature(X = "Expression", Y = "ANY"), function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  if(FUN != "*" || make.dimnames)
    stop("Unimplemented")
  Kron(X, Y)
})

#' @docType methods
#' @rdname kronecker
#' @export
setMethod("kronecker", signature(X = "ANY", Y = "Expression"), function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  if(FUN != "*" || make.dimnames)
    stop("Unimplemented")
  Kron(X, Y)
})

#' @docType methods
#' @rdname kronecker
#' @export
setMethod("%x%", signature(X = "Expression", Y = "ANY"), function(X, Y) { Kron(lh_exp = X, rh_exp = Y) })

#' @docType methods
#' @rdname kronecker
#' @export
setMethod("%x%", signature(X = "ANY", Y = "Expression"), function(X, Y) { Kron(lh_exp = X, rh_exp = Y) })
