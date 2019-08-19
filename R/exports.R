# =========================
#     Scalar functions
# =========================
#'
#' Unity Resolvent
#' 
#' The unity resolvent of a positive matrix. For an elementwise positive matrix \eqn{X}, this atom represents \eqn{(I - X)^{-1}},
#' and it enforces the constraint that the spectral radius of \eqn{X} is at most 1.
#' 
#' This atom is log-log convex.
#' 
#' @param X An \linkS4class{Expression} or positive square matrix.
#' @return An \linkS4class{Expression} representing the unity resolvent of the input.
#' @docType methods
#' @name eye_minus_inv
#' @rdname eye_minus_inv
#' @export
eye_minus_inv <- EyeMinusInv

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
#' @param max_denom (Optional) The maximum denominator to use in approximating \code{p/sum(p)} with \code{w}. If \code{w} is not an exact representation, increasing \code{max_denom} may offer a more accurate representation, at the cost of requiring more convex inequalities to represent the geometric mean. Defaults to 1024.
#' @return An \linkS4class{Expression} representing the geometric mean of the input.
#' @examples
#' x <- Variable(2)
#' cost <- geo_mean(x)
#' prob <- Problem(Maximize(cost), list(sum(x) <= 1))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#'
#' x <- Variable(5)
#' p <- c(0.07, 0.12, 0.23, 0.19, 0.39)
#' prob <- Problem(Maximize(geo_mean(x,p)), list(p_norm(x) <= 1))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' @docType methods
#' @name geo_mean
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
#' x <- Variable()
#' prob <- Problem(Maximize(harmonic_mean(x)), list(x >= 0, x <= 5))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' @docType methods
#' @name harmonic_mean
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
#' A <- Variable(2,2)
#' prob <- Problem(Minimize(lambda_max(A)), list(A >= 2))
#' result <- solve(prob)
#' result$value
#' result$getValue(A)
#'
#' obj <- Maximize(A[2,1] - A[1,2])
#' prob <- Problem(obj, list(lambda_max(A) <= 100, A[1,1] == 2, A[2,2] == 2, A[2,1] == 2))
#' result <- solve(prob)
#' result$value
#' result$getValue(A)
#' @docType methods
#' @name lambda_max
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
#' A <- Variable(2,2)
#' val <- cbind(c(5,7), c(7,-3))
#' prob <- Problem(Maximize(lambda_min(A)), list(A == val))
#' result <- solve(prob)
#' result$value
#' result$getValue(A)
#' @docType methods
#' @name lambda_min
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
#' @examples
#' C <- Variable(3,3)
#' val <- cbind(c(1,2,3), c(2,4,5), c(3,5,6))
#' prob <- Problem(Minimize(lambda_sum_largest(C,2)), list(C == val))
#' result <- solve(prob)
#' result$value
#' result$getValue(C)
#' @docType methods
#' @name lambda_sum_largest
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
#' @examples
#' C <- Variable(3,3)
#' val <- cbind(c(1,2,3), c(2,4,5), c(3,5,6))
#' prob <- Problem(Maximize(lambda_sum_smallest(C,2)), list(C == val))
#' result <- solve(prob)
#' result$value
#' result$getValue(C)
#' @docType methods
#' @name lambda_sum_smallest
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
#' @examples
#' x <- t(data.frame(c(0.55, 0.25, -0.2, -0.25, -0.0, 0.4),
#'                   c(0.0, 0.35, 0.2, -0.1, -0.3, -0.2)))
#' n <- nrow(x)
#' m <- ncol(x)
#'
#' A <- Variable(n,n)
#' b <- Variable(n)
#' obj <- Maximize(log_det(A))
#' constr <- lapply(1:m, function(i) { p_norm(A %*% as.matrix(x[,i]) + b) <= 1 })
#' prob <- Problem(obj, constr)
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @name log_det
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
#' @examples
#' A <- Variable(2,2)
#' val <- cbind(c(5,7), c(0,-3))
#' prob <- Problem(Minimize(log_sum_exp(A)), list(A == val))
#' result <- solve(prob)
#' result$getValue(A)
#' @docType methods
#' @name log_sum_exp
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
#' @examples
#' \dontrun{
#' m <- 100
#' n <- 80
#' r <- 70
#'
#' A <- matrix(stats::rnorm(m*n), nrow = m, ncol = n)
#' b <- matrix(stats::rnorm(m), nrow = m, ncol = 1)
#' G <- matrix(stats::rnorm(r*n), nrow = r, ncol = n)
#' h <- matrix(stats::rnorm(r), nrow = r, ncol = 1)
#'
#' # ||Ax-b||^2 = x^T (A^T A) x - 2(A^T b)^T x + ||b||^2
#' P <- t(A) %*% A
#' q <- -2 * t(A) %*% b
#' r <- t(b) %*% b
#' Pinv <- base::solve(P)
#'
#' x <- Variable(n)
#' obj <- matrix_frac(x, Pinv) + t(q) %*% x + r
#' constr <- list(G %*% x == h)
#' prob <- Problem(Minimize(obj), constr)
#' result <- solve(prob)
#' result$value
#' }
#' @docType methods
#' @name matrix_frac
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
#' @examples
#' x <- Variable(2)
#' val <- matrix(c(-5,-10))
#' prob <- Problem(Minimize(max_entries(x)), list(x == val))
#' result <- solve(prob)
#' result$value
#'
#' A <- Variable(2,2)
#' val <- rbind(c(-5,2), c(-3,1))
#' prob <- Problem(Minimize(max_entries(A, axis = 1)[2,1]), list(A == val))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @name max_entries
#' @aliases max
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
#' @examples
#' A <- Variable(2,2)
#' val <- cbind(c(-5,2), c(-3,1))
#' prob <- Problem(Maximize(min_entries(A)), list(A == val))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @name min_entries
#' @aliases min
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
#' @return An \linkS4class{Expression} representing the \eqn{l_{p,q}} norm of the input.
#' @examples
#' A <- Variable(2,2)
#' val <- cbind(c(3,3), c(4,4))
#' prob <- Problem(Minimize(mixed_norm(A,2,1)), list(A == val))
#' result <- solve(prob)
#' result$value
#' result$getValue(A)
#'
#' val <- cbind(c(1,4), c(5,6))
#' prob <- Problem(Minimize(mixed_norm(A,1,Inf)), list(A == val))
#' result <- solve(prob)
#' result$value
#' result$getValue(A)
#' @docType methods
#' @name mixed_norm
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
#' a <- Variable()
#' prob <- Problem(Minimize(norm1(a)), list(a <= -2))
#' result <- solve(prob)
#' result$value
#' result$getValue(a)
#'
#' prob <- Problem(Maximize(-norm1(a)), list(a <= -2))
#' result <- solve(prob)
#' result$value
#' result$getValue(a)
#'
#' x <- Variable(2)
#' z <- Variable(2)
#' prob <- Problem(Minimize(norm1(x - z) + 5), list(x >= c(2,3), z <= c(-1,-4)))
#' result <- solve(prob)
#' result$value
#' result$getValue(x[1] - z[1])
#' @docType methods
#' @name norm1
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
#' a <- Variable()
#' prob <- Problem(Minimize(norm2(a)), list(a <= -2))
#' result <- solve(prob)
#' result$value
#' result$getValue(a)
#'
#' prob <- Problem(Maximize(-norm2(a)), list(a <= -2))
#' result <- solve(prob)
#' result$value
#' result$getValue(a)
#'
#' x <- Variable(2)
#' z <- Variable(2)
#' prob <- Problem(Minimize(norm2(x - z) + 5), list(x >= c(2,3), z <= c(-1,-4)))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' result$getValue(z)
#'
#' prob <- Problem(Minimize(norm2(t(x - z)) + 5), list(x >= c(2,3), z <= c(-1,-4)))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' result$getValue(z)
#' @docType methods
#' @name norm2
#' @rdname norm2
#' @export
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
#' a <- Variable()
#' b <- Variable()
#' c <- Variable()
#'
#' prob <- Problem(Minimize(norm_inf(a)), list(a >= 2))
#' result <- solve(prob)
#' result$value
#' result$getValue(a)
#'
#' prob <- Problem(Minimize(3*norm_inf(a + 2*b) + c), list(a >= 2, b <= -1, c == 3))
#' result <- solve(prob)
#' result$value
#' result$getValue(a + 2*b)
#' result$getValue(c)
#'
#' prob <- Problem(Maximize(-norm_inf(a)), list(a <= -2))
#' result <- solve(prob)
#' result$value
#' result$getValue(a)
#'
#' x <- Variable(2)
#' z <- Variable(2)
#' prob <- Problem(Minimize(norm_inf(x - z) + 5), list(x >= c(2,3), z <= c(-1,-4)))
#' result <- solve(prob)
#' result$value
#' result$getValue(x[1] - z[1])
#' @docType methods
#' @name norm_inf
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
#' @examples
#' C <- Variable(3,3)
#' val <- cbind(3:5, 6:8, 9:11)
#' prob <- Problem(Minimize(norm_nuc(C)), list(C == val))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @name norm_nuc
#' @rdname norm_nuc
#' @export
norm_nuc <- NormNuc

#'
#' Difference on Restricted Domain
#' 
#' The difference \eqn{1 - x} with domain \eqn{\{x : 0 < x < 1\}}.
#'
#' This atom is log-log concave.
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @return An \linkS4class{Expression} representing one minus the input restricted to \eqn{(0,1)}.
#' @docType methods
#' @name one_minus_pos
#' @rdname one_minus_pos
#' @export
one_minus_pos <- OneMinusPos

#'
#' Perron-Frobenius Eigenvalue
#' 
#' The Perron-Frobenius eigenvalue of a positive matrix.
#' 
#' For an elementwise positive matrix \eqn{X}, this atom represents its spectral radius, i.e., the magnitude of its largest eigenvalue.
#' Because \eqn{X} is positive, the spectral radius equals its largest eigenvalue, which is guaranteed to be positive.
#' 
#' This atom is log-log convex.
#' @param X An \linkS4class{Expression} or positive square matrix.
#' @return An \linkS4class{Expression} representing the largest eigenvalue of the input.
#' @docType methods
#' @name pf_eigenvalue
#' @rdname pf_eigenvalue
#' @export
pf_eigenvalue <- PfEigenvalue

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
#' @param max_denom The maximum denominator considered in forming a rational approximation for \eqn{p}.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \linkS4class{Expression} representing the p-norm of the input.
#' @examples
#' x <- Variable(3)
#' prob <- Problem(Minimize(p_norm(x,2)))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#'
#' prob <- Problem(Minimize(p_norm(x,Inf)))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#'
#' a <- c(1.0, 2, 3)
#' prob <- Problem(Minimize(p_norm(x,1.6)), list(t(x) %*% a >= 1))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#'
#' prob <- Problem(Minimize(sum(abs(x - a))), list(p_norm(x,-1) >= 0))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' @docType methods
#' @name p_norm
#' @rdname p_norm
#' @export
p_norm <- Pnorm

#'
#' Product of Entries
#'
#' The product of entries in a vector or matrix.
#' 
#' This atom is log-log affine, but it is neither convex nor concave.
#'
#' @param expr An \linkS4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \linkS4class{Expression} representing the product of the entries of the input.
#' @docType methods
#' @name prod_entries
#' @aliases prod
#' @rdname prod_entries
#' @export
prod_entries <- ProdEntries

#'
#' Quadratic Form
#'
#' The quadratic form, \eqn{x^TPx}.
#'
#' @param x An \linkS4class{Expression} or vector.
#' @param P An \linkS4class{Expression} or matrix.
#' @return An \linkS4class{Expression} representing the quadratic form evaluated at the input.
#' @examples
#' x <- Variable(2)
#' P <- rbind(c(4,0), c(0,9))
#' prob <- Problem(Minimize(quad_form(x,P)), list(x >= 1))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#'
#' A <- Variable(2,2)
#' c <- c(1,2)
#' prob <- Problem(Minimize(quad_form(c,A)), list(A >= 1))
#' result <- solve(prob)
#' result$value
#' result$getValue(A)
#' @docType methods
#' @name quad_form
#' @rdname quad_form
#' @export
quad_form <- function(x, P) {
  # x^T P x
  x <- as.Constant(x)
  P <- as.Constant(P)
  
  # Check dimensions.
  P_dim <- dim(P)
  if(ndim(P) != 2 || P_dim[1] != P_dim[2] || max(nrow(x), 1) != P_dim[1])
    stop("Invalid dimensions for arguments.")
  
  # P cannot be a parameter.
  if(is_constant(x))
    Conj(t(x)) %*% P %*% x
  else if(is_constant(P))
    QuadForm(x, P)
  else
    stop("At least one argument to QuadForm must be constant.")
}

#'
#' Quadratic over Linear
#'
#' \eqn{\sum_{i,j} X_{i,j}^2/y}.
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param y A scalar \linkS4class{Expression} or numeric constant.
#' @return An \linkS4class{Expression} representing the quadratic over linear function value evaluated at the input.
#' @examples
#' x <- Variable(3,2)
#' y <- Variable()
#' val <- cbind(c(-1,2,-2), c(-1,2,-2))
#' prob <- Problem(Minimize(quad_over_lin(x,y)), list(x == val, y <= 2))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' result$getValue(y)
#' @docType methods
#' @name quad_over_lin
#' @rdname quad_over_lin
#' @export
quad_over_lin <- QuadOverLin

#'
#' Sum of Entries
#'
#' The sum of entries in a vector or matrix.
#'
#' @param expr An \linkS4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \linkS4class{Expression} representing the sum of the entries of the input.
#' @examples
#' x <- Variable(2)
#' prob <- Problem(Minimize(sum_entries(x)), list(t(x) >= matrix(c(1,2), nrow = 1, ncol = 2)))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#'
#' C <- Variable(3,2)
#' prob <- Problem(Maximize(sum_entries(C)), list(C[2:3,] <= 2, C[1,] == 1))
#' result <- solve(prob)
#' result$value
#' result$getValue(C)
#' @docType methods
#' @name sum_entries
#' @aliases sum
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
#' @examples
#' m <- 300
#' n <- 9
#' X <- matrix(stats::rnorm(m*n), nrow = m, ncol = n)
#' X <- cbind(rep(1,m), X)
#' b <- c(0, 0.8, 0, 1, 0.2, 0, 0.4, 1, 0, 0.7)
#' y <- X %*% b + stats::rnorm(m)
#'
#' beta <- Variable(n+1)
#' obj <- sum_largest((y - X %*% beta)^2, 100)
#' prob <- Problem(Minimize(obj))
#' result <- solve(prob)
#' result$getValue(beta)
#' @docType methods
#' @name sum_largest
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
#' @examples
#' m <- 300
#' n <- 9
#' X <- matrix(stats::rnorm(m*n), nrow = m, ncol = n)
#' X <- cbind(rep(1,m), X)
#' b <- c(0, 0.8, 0, 1, 0.2, 0, 0.4, 1, 0, 0.7)
#' factor <- 2*rbinom(m, size = 1, prob = 0.8) - 1
#' y <- factor * (X %*% b) + stats::rnorm(m)
#'
#' beta <- Variable(n+1)
#' obj <- sum_smallest(y - X %*% beta, 200)
#' prob <- Problem(Maximize(obj), list(0 <= beta, beta <= 1))
#' result <- solve(prob)
#' result$getValue(beta)
#' @docType methods
#' @name sum_smallest
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
#' @examples
#' m <- 30
#' n <- 20
#' A <- matrix(stats::rnorm(m*n), nrow = m, ncol = n)
#' b <- matrix(stats::rnorm(m), nrow = m, ncol = 1)
#'
#' x <- Variable(n)
#' obj <- Minimize(sum_squares(A %*% x - b))
#' constr <- list(0 <= x, x <= 1)
#' prob <- Problem(obj, constr)
#' result <- solve(prob)
#'
#' result$value
#' result$getValue(x)
#' result$getDualValue(constr[[1]])
#' @docType methods
#' @name sum_squares
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
#' @examples
#' C <- Variable(3,3)
#' val <- cbind(3:5, 6:8, 9:11)
#' prob <- Problem(Maximize(matrix_trace(C)), list(C == val))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @name matrix_trace
#' @aliases trace tr
#' @rdname matrix_trace
#' @export
matrix_trace <- Trace

#'
#' Total Variation
#'
#' The total variation of a vector, matrix, or list of matrices. Uses L1 norm of discrete gradients for vectors and L2 norm of discrete gradients for matrices.
#'
#' @param value An \linkS4class{Expression}, vector, or matrix.
#' @param ... (Optional) \linkS4class{Expression} objects or numeric constants that extend the third dimension of value.
#' @return An \linkS4class{Expression} representing the total variation of the input.
#' @examples
#' rows <- 10
#' cols <- 10
#' Uorig <- matrix(sample(0:255, size = rows * cols, replace = TRUE), nrow = rows, ncol = cols)
#'
#' # Known is 1 if the pixel is known, 0 if the pixel was corrupted
#' Known <- matrix(0, nrow = rows, ncol = cols)
#' for(i in 1:rows) {
#'    for(j in 1:cols) {
#'       if(stats::runif(1) > 0.7)
#'          Known[i,j] <- 1
#'    }
#' }
#' Ucorr <- Known %*% Uorig
#'
#' # Recover the original image using total variation in-painting
#' U <- Variable(rows, cols)
#' obj <- Minimize(tv(U))
#' constraints <- list(Known * U == Known * Ucorr)
#' prob <- Problem(obj, constraints)
#' result <- solve(prob, solver = "SCS")
#' result$getValue(U)
#' @docType methods
#' @name tv
#' @aliases total_variation
#' @rdname tv
#' @export
tv <- TotalVariation

#' @param ... Numeric scalar, vector, matrix, or \linkS4class{Expression} objects.
#' @param na.rm (Unimplemented) A logical value indicating whether missing values should be removed.
#' @docType methods
#' @rdname max_entries
#' @method max Expression
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

#' @param ... Numeric scalar, vector, matrix, or \linkS4class{Expression} objects.
#' @param na.rm (Unimplemented) A logical value indicating whether missing values should be removed.
#' @docType methods
#' @rdname min_entries
#' @method min Expression
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
#' @param x An \linkS4class{Expression}.
#' @param type A character indicating the type of norm desired.
#' \itemize{
#'    \item "O", "o" or "1" specifies the 1-norm (maximum absolute column sum).
#'    \item "I" or "i" specifies the infinity-norm (maximum absolute row sum).
#'    \item "F" or "f" specifies the Frobenius norm (Euclidean norm of the vectorized \code{x}).
#'    \item "M" or "m" specifies the maximum modulus of all the elements in \code{x}.
#'    \item "2" specifies the spectral norm, which is the largest singular value of \code{x}.
#' }
#' @return An \linkS4class{Expression} representing the norm of the input.
#' @seealso The \code{\link{p_norm}} function calculates the vector p-norm.
#' @examples
#' C <- Variable(3,2)
#' val <- Constant(rbind(c(1,2), c(3,4), c(5,6)))
#' prob <- Problem(Minimize(norm(C, "F")), list(C == val))
#' result <- solve(prob, solver = "SCS")
#' result$value
#' @docType methods
#' @aliases norm
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

#'
#' Matrix Norm (Alternative)
#'
#' A wrapper on the different norm atoms. This is different from the standard "norm" method in the R base package.
#' If \code{p = 2}, \code{axis = NA}, and \code{x} is a matrix, this returns the maximium singular value.
#'
#' @param x An \linkS4class{Expression} or numeric constant representing a vector or matrix.
#' @param p The type of norm. May be a number (p-norm), "inf" (infinity-norm), "nuc" (nuclear norm), or "fro" (Frobenius norm). The default is \code{p = 2}.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \linkS4class{Expression} representing the norm.
#' @seealso \link[CVXR]{norm}
#' @docType methods
#' @name cvxr_norm
#' @rdname cvxr_norm
#' @export
cvxr_norm <- Norm

#' @param ... Numeric scalar, vector, matrix, or \linkS4class{Expression} objects.
#' @param na.rm (Unimplemented) A logical value indicating whether missing values should be removed.
#' @docType methods
#' @rdname sum_entries
#' @method sum Expression
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

#' @param ... Numeric scalar, vector, matrix, or \linkS4class{Expression} objects.
#' @param na.rm (Unimplemented) A logical value indicating whether missing values should be removed.
#' @docType methods
#' @rdname prod_entries
#' @method prod Expression
#' @export
prod.Expression <- function(..., na.rm = FALSE) {
  if(na.rm)
    warning("na.rm is unimplemented for Expression objects")
  
  vals <- list(...)
  is_expr <- sapply(vals, function(v) { is(v, "Expression") })
  sum_expr <- lapply(vals[is_expr], function(expr) { ProdEntries(expr = expr) })
  if(all(is_expr))
    Reduce("*", sum_expr)
  else {
    sum_num <- sum(sapply(vals[!is_expr], function(v) { prod(v, na.rm = na.rm) }))
    Reduce("*", sum_expr) + sum_num
  }
}

#'
#' Arithmetic Mean
#'
#' The arithmetic mean of an expression.
#'
#' @param x An \linkS4class{Expression} object.
#' @param trim (Unimplemented) The fraction (0 to 0.5) of observations to be trimmed from each end of \eqn{x} before the mean is computed.
#' @param na.rm (Unimplemented) A logical value indicating whether missing values should be removed.
#' @param ... (Unimplemented) Optional arguments.
#' @return An \linkS4class{Expression} representing the mean of the input.
#' @examples
#' A <- Variable(2,2)
#' val <- cbind(c(-5,2), c(-3,1))
#' prob <- Problem(Minimize(mean(A)), list(A == val))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @aliases mean
#' @rdname mean
#' @method mean Expression
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
#' @examples
#' x <- Variable(5)
#' obj <- Maximize(sum(entr(x)))
#' prob <- Problem(obj, list(sum(x) == 1))
#' result <- solve(prob)
#' result$getValue(x)
#' @docType methods
#' @name entr
#' @aliases entropy
#' @rdname entr
#' @export
entr <- Entr

#'
#' Huber Function
#'
#' The elementwise Huber function, \eqn{Huber(x, M) = }
#' \itemize{
#'   \item{\eqn{2M|x|-M^2}}{for \eqn{|x| \geq |M|}}
#'    \item{\eqn{|x|^2}}{for \eqn{|x| \leq |M|.}}
#'  }
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param M (Optional) A positive scalar value representing the threshold. Defaults to 1.
#' @return An \linkS4class{Expression} representing the Huber function evaluated at the input.
#' @examples
#' n <- 10
#' m <- 450
#' p <- 0.1    # Fraction of responses with sign flipped
#'
#' # Generate problem data
#' beta_true <- 5*matrix(stats::rnorm(n), nrow = n)
#' X <- matrix(stats::rnorm(m*n), nrow = m, ncol = n)
#' y_true <- X %*% beta_true
#' eps <- matrix(stats::rnorm(m), nrow = m)
#'
#' # Randomly flip sign of some responses
#' factor <- 2*rbinom(m, size = 1, prob = 1-p) - 1
#' y <- factor * y_true + eps
#'
#' # Huber regression
#' beta <- Variable(n)
#' obj <- sum(huber(y - X %*% beta, 1))
#' prob <- Problem(Minimize(obj))
#' result <- solve(prob)
#' result$getValue(beta)
#' @docType methods
#' @name huber
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
#' @examples
#' A <- Variable(2,2)
#' val <- cbind(c(1,2), c(3,4))
#' prob <- Problem(Minimize(inv_pos(A)[1,2]), list(A == val))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @name inv_pos
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
#' @examples
#' n <- 5
#' alpha <- seq(10, n-1+10)/n
#' beta <- seq(10, n-1+10)/n
#' P_tot <- 0.5
#' W_tot <- 1.0
#'
#' P <- Variable(n)
#' W <- Variable(n)
#' R <- kl_div(alpha*W, alpha*(W + beta*P)) - alpha*beta*P
#' obj <- sum(R)
#' constr <- list(P >= 0, W >= 0, sum(P) == P_tot, sum(W) == W_tot)
#' prob <- Problem(Minimize(obj), constr)
#' result <- solve(prob)
#'
#' result$value
#' result$getValue(P)
#' result$getValue(W)
#' @docType methods
#' @name kl_div
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
#' @examples
#' n <- 20
#' m <- 1000
#' sigma <- 45
#'
#' beta_true <- stats::rnorm(n)
#' idxs <- sample(n, size = 0.8*n, replace = FALSE)
#' beta_true[idxs] <- 0
#' X <- matrix(stats::rnorm(m*n, 0, 5), nrow = m, ncol = n)
#' y <- sign(X %*% beta_true + stats::rnorm(m, 0, sigma))
#'
#' beta <- Variable(n)
#' X_sign <- apply(X, 2, function(x) { ifelse(y <= 0, -1, 1) * x })
#' obj <- -sum(logistic(-X[y <= 0,] %*% beta)) - sum(logistic(X[y == 1,] %*% beta))
#' prob <- Problem(Maximize(obj))
#' result <- solve(prob)
#'
#' log_odds <- result$getValue(X %*% beta)
#' beta_res <- result$getValue(beta)
#' y_probs <- 1/(1 + exp(-X %*% beta_res))
#' log(y_probs/(1 - y_probs))
#' @docType methods
#' @name logistic
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
#' @examples
#' c <- matrix(c(1,-1))
#' prob <- Problem(Minimize(max_elemwise(t(c), 2, 2 + t(c))[2]))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @name max_elemwise
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
#' @examples
#' a <- cbind(c(-5,2), c(-3,-1))
#' b <- cbind(c(5,4), c(-1,2))
#' prob <- Problem(Minimize(min_elemwise(a, 0, b)[1,2]))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @name min_elemwise
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
#' A <- Variable(2,2)
#' c <- cbind(c(1,-1), c(2,-2))
#' expr <- multiply(c, A)
#' obj <- Minimize(norm_inf(expr))
#' prob <- Problem(obj, list(A == 5))
#' result <- solve(prob)
#' result$value
#' result$getValue(expr)
#' @docType methods
#' @name multiply
#' @aliases *
#' @rdname multiply
#' @export
multiply <- Multiply

#'
#' Elementwise Negative
#'
#' The elementwise absolute negative portion of an expression, \eqn{-\min(x_i,0)}. This is equivalent to \code{-min_elemwise(x,0)}.
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @return An \linkS4class{Expression} representing the negative portion of the input.
#' @examples
#' x <- Variable(2)
#' val <- matrix(c(-3,3))
#' prob <- Problem(Minimize(neg(x)[1]), list(x == val))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @name neg
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
#' @examples
#' x <- Variable(2)
#' val <- matrix(c(-3,2))
#' prob <- Problem(Minimize(pos(x)[1]), list(x == val))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @name pos
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
#' For \eqn{p < 0} and \eqn{f(x) = }
#' \itemize{
#'   \item{\eqn{x^p}}{ for \eqn{x > 0}}
#'   \item{\eqn{+\infty}}{\eqn{x \leq 0}}
#' }, this function is convex, decreasing, and positive.
#' For \eqn{0 < p < 1} and \eqn{f(x) =}
#' \itemize{
#'   \item{\eqn{x^p}}{ for \eqn{x \geq 0}}
#'   \item{\eqn{-\infty}}{\eqn{x < 0}}
#' }, this function is concave, increasing, and positivea.
#' For \eqn{p > 1, p \neq 2,4,8,\ldots} and \eqn{f(x) = }
#' \itemize{
#'   \item{\eqn{x^p}}{ for \eqn{x \geq 0}}
#'   \item{\eqn{+\infty}}{\eqn{x < 0}}
#' }, this function is convex, increasing, and positive.
#'
#' @param x An \linkS4class{Expression}, vector, or matrix.
#' @param p A scalar value indicating the exponential power.
#' @param max_denom The maximum denominator considered in forming a rational approximation of \code{p}.
#' @examples
#' \dontrun{
#' x <- Variable()
#' prob <- Problem(Minimize(power(x,1.7) + power(x,-2.3) - power(x,0.45)))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' }
#' @docType methods
#' @name power
#' @aliases ^
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
#' @examples
#' \dontrun{
#' A <- Variable(2,2)
#' val <- cbind(c(-5,2), c(-3,1))
#' prob <- Problem(Minimize(scalene(A,2,3)[1,1]), list(A == val))
#' result <- solve(prob)
#' result$value
#' result$getValue(scalene(A, 0.7, 0.3))
#' }
#' @docType methods
#' @name scalene
#' @rdname scalene
#' @export
scalene <- Scalene

#'
#' Absolute Value
#'
#' The elementwise absolute value.
#'
#' @param x An \linkS4class{Expression}.
#' @return An \linkS4class{Expression} representing the absolute value of the input.
#' @examples
#' A <- Variable(2,2)
#' prob <- Problem(Minimize(sum(abs(A))), list(A <= -2))
#' result <- solve(prob)
#' result$value
#' result$getValue(A)
#' @docType methods
#' @aliases abs
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
#' @examples
#' x <- Variable(5)
#' obj <- Minimize(sum(exp(x)))
#' prob <- Problem(obj, list(sum(x) == 1))
#' result <- solve(prob)
#' result$getValue(x)
#' @docType methods
#' @aliases exp
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
#' @examples
#' # Log in objective
#' x <- Variable(2)
#' obj <- Maximize(sum(log(x)))
#' constr <- list(x <= matrix(c(1, exp(1))))
#' prob <- Problem(obj, constr)
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#'
#' # Log in constraint
#' obj <- Minimize(sum(x))
#' constr <- list(log2(x) >= 0, x <= matrix(c(1,1)))
#' prob <- Problem(obj, constr)
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#'
#' # Index into log
#' obj <- Maximize(log10(x)[2])
#' constr <- list(x <= matrix(c(1, exp(1))))
#' prob <- Problem(obj, constr)
#' result <- solve(prob)
#' result$value
#'
#' # Scalar log
#' obj <- Maximize(log1p(x[2]))
#' constr <- list(x <= matrix(c(1, exp(1))))
#' prob <- Problem(obj, constr)
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @aliases log log10 log2 log1p
#' @rdname log
#' @export
setMethod("log", "Expression", function(x, base = base::exp(1)) {
  if(base == base::exp(1))
    Log(x = x)
  else
    Log(x = x)/base::log(base)
})

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
#' A <- Variable(2,2)
#' val <- cbind(c(2,4), c(16,1))
#' prob <- Problem(Maximize(sqrt(A)[1,2]), list(A == val))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @aliases sqrt
#' @rdname sqrt
#' @export
setMethod("sqrt", "Expression", function(x) { Power(x = x, p = 0.5) })

#'
#' Square
#'
#' The elementwise square.
#'
#' @param x An \linkS4class{Expression}.
#' @return An \linkS4class{Expression} representing the square of the input.
#' A <- Variable(2,2)
#' val <- cbind(c(2,4), c(16,1))
#' prob <- Problem(Minimize(square(A)[1,2]), list(A == val))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @aliases square
#' @rdname square
#' @export
setMethod("square", "Expression", function(x) { Power(x = x, p = 2) })

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
#' @examples
#' x <- Variable()
#' expr <- bmat(list(list(matrix(1, nrow = 3, ncol = 1), matrix(2, nrow = 3, ncol = 2)),
#'                 list(matrix(3, nrow = 1, ncol = 2), x)
#'              ))
#' prob <- Problem(Minimize(sum_entries(expr)), list(x >= 0))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @name bmat
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
#' h <- matrix(stats::rnorm(2), nrow = 2, ncol = 1)
#' prob <- Problem(Minimize(sum(conv(h, x))))
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' @docType methods
#' @name conv
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
#' x <- Variable(2)
#' y <- Variable(3)
#' c <- matrix(1, nrow = 1, ncol = 5)
#' prob <- Problem(Minimize(c %*% t(hstack(t(x), t(y)))), list(x == c(1,2), y == c(3,4,5)))
#' result <- solve(prob)
#' result$value
#'
#' c <- matrix(1, nrow = 1, ncol = 4)
#' prob <- Problem(Minimize(c %*% t(hstack(t(x), t(x)))), list(x == c(1,2)))
#' result <- solve(prob)
#' result$value
#'
#' A <- Variable(2,2)
#' C <- Variable(3,2)
#' c <- matrix(1, nrow = 2, ncol = 2)
#' prob <- Problem(Minimize(sum_entries(hstack(t(A), t(C)))), list(A >= 2*c, C == -2))
#' result <- solve(prob)
#' result$value
#' result$getValue(A)
#'
#' D <- Variable(3,3)
#' expr <- hstack(C, D)
#' obj <- expr[1,2] + sum(hstack(expr, expr))
#' constr <- list(C >= 0, D >= 0, D[1,1] == 2, C[1,2] == 3)
#' prob <- Problem(Minimize(obj), constr)
#' result <- solve(prob)
#' result$value
#' result$getValue(C)
#' result$getValue(D)
#' @docType methods
#' @name hstack
#' @rdname hstack
#' @export
hstack <- HStack

#'
#' Reshape an Expression
#'
#' This function vectorizes an expression, then unvectorizes it into a new shape. Entries are stored in column-major order.
#'
#' @param expr An \linkS4class{Expression}, vector, or matrix.
#' @param new_dim The new dimensions.
#' @return An \linkS4class{Expression} representing the reshaped input.
#' @examples
#' x <- Variable(4)
#' mat <- cbind(c(1,-1), c(2,-2))
#' vec <- matrix(1:4)
#' expr <- reshape_expr(x,c(2,2))
#' obj <- Minimize(sum(mat %*% expr))
#' prob <- Problem(obj, list(x == vec))
#' result <- solve(prob)
#' result$value
#'
#' A <- Variable(2,2)
#' c <- 1:4
#' expr <- reshape_expr(A,c(4,1))
#' obj <- Minimize(t(expr) %*% c)
#' constraints <- list(A == cbind(c(-1,-2), c(3,4)))
#' prob <- Problem(obj, constraints)
#' result <- solve(prob)
#' result$value
#' result$getValue(expr)
#' result$getValue(reshape_expr(expr,c(2,2)))
#'
#' C <- Variable(3,2)
#' expr <- reshape_expr(C,c(2,3))
#' mat <- rbind(c(1,-1), c(2,-2))
#' C_mat <- rbind(c(1,4), c(2,5), c(3,6))
#' obj <- Minimize(sum(mat %*% expr))
#' prob <- Problem(obj, list(C == C_mat))
#' result <- solve(prob)
#' result$value
#' result$getValue(expr)
#'
#' a <- Variable()
#' c <- cbind(c(1,-1), c(2,-2))
#' expr <- reshape_expr(c * a,c(1,4))
#' obj <- Minimize(expr %*% (1:4))
#' prob <- Problem(obj, list(a == 2))
#' result <- solve(prob)
#' result$value
#' result$getValue(expr)
#'
#' expr <- reshape_expr(c * a,c(4,1))
#' obj <- Minimize(t(expr) %*% (1:4))
#' prob <- Problem(obj, list(a == 2))
#' result <- solve(prob)
#' result$value
#' result$getValue(expr)
#' @docType methods
#' @name reshape_expr
#' @aliases reshape
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
#' @examples
#' C <- Variable(3,2)
#' val <- rbind(c(1,2), c(3,4), c(5,6))
#' obj <- sigma_max(C)
#' constr <- list(C == val)
#' prob <- Problem(Minimize(obj), constr)
#' result <- solve(prob, solver = "SCS")
#' result$value
#' result$getValue(C)
#' @docType methods
#' @name sigma_max
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
#' @examples
#' C <- Variable(3,3)
#' val <- cbind(3:5, 6:8, 9:11)
#' prob <- Problem(Maximize(upper_tri(C)[3,1]), list(C == val))
#' result <- solve(prob)
#' result$value
#' result$getValue(upper_tri(C))
#' @docType methods
#' @name upper_tri
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
#' A <- Variable(2,2)
#' c <- 1:4
#' expr <- vec(A)
#' obj <- Minimize(t(expr) %*% c)
#' constraints <- list(A == cbind(c(-1,-2), c(3,4)))
#' prob <- Problem(obj, constraints)
#' result <- solve(prob)
#' result$value
#' result$getValue(expr)
#' @docType methods
#' @name vec
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
#' x <- Variable(2)
#' y <- Variable(3)
#' c <- matrix(1, nrow = 1, ncol = 5)
#' prob <- Problem(Minimize(c %*% vstack(x, y)), list(x == c(1,2), y == c(3,4,5)))
#' result <- solve(prob)
#' result$value
#'
#' c <- matrix(1, nrow = 1, ncol = 4)
#' prob <- Problem(Minimize(c %*% vstack(x, x)), list(x == c(1,2)))
#' result <- solve(prob)
#' result$value
#'
#' A <- Variable(2,2)
#' C <- Variable(3,2)
#' c <- matrix(1, nrow = 2, ncol = 2)
#' prob <- Problem(Minimize(sum(vstack(A, C))), list(A >= 2*c, C == -2))
#' result <- solve(prob)
#' result$value
#'
#' B <- Variable(2,2)
#' c <- matrix(1, nrow = 1, ncol = 2)
#' prob <- Problem(Minimize(sum(vstack(c %*% A, c %*% B))), list(A >= 2, B == -2))
#' result <- solve(prob)
#' result$value
#' @docType methods
#' @name vstack
#' @rdname vstack
#' @export
vstack <- VStack

#'
#' Cumulative Sum
#'
#' The cumulative sum, \eqn{\sum_{i=1}^k x_i} for \eqn{k=1,\ldots,n}.
#' When calling \code{cumsum}, matrices are automatically flattened into column-major order before the sum is taken.
#'
#' @param x,expr An \linkS4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, and \code{2} indicates columns. The default is \code{2}.
#' @examples
#' val <- cbind(c(1,2), c(3,4))
#' value(cumsum(Constant(val)))
#' value(cumsum_axis(Constant(val)))
#'
#' x <- Variable(2,2)
#' prob <- Problem(Minimize(cumsum(x)[4]), list(x == val))
#' result <- solve(prob)
#' result$value
#' result$getValue(cumsum(x))
#' @docType methods
#' @name cumsum_axis
#' @aliases cumsum_axis cumsum
#' @rdname cumsum_axis
#' @export
cumsum_axis <- CumSum

#' @docType methods
#' @rdname cumsum_axis
#' @export
setMethod("cumsum", signature(x = "Expression"), function(x) { CumSum(expr = Vec(x)) })

#'
#' Matrix Diagonal
#'
#' Extracts the diagonal from a matrix or makes a vector into a diagonal matrix.
#'
#' @param x An \linkS4class{Expression}, vector, or square matrix.
#' @param nrow,ncol (Optional) Dimensions for the result when \code{x} is not a matrix.
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
#' @aliases diag
#' @rdname diag
#' @export
setMethod("diag", signature(x = "Expression"), function(x = 1, nrow, ncol) {
    missing_nrow <- missing(nrow)
    missing_ncol <- missing(ncol)
    if (missing_nrow && missing_ncol) {
        Diag(x)
    } else if (is_matrix(x) & ((!missing_nrow) || (!missing_ncol))) {
        stop("'nrow' or 'ncol' cannot be specified when 'x' is a matrix")
    } else {
        expr <- as.Constant(x)
        n <- length(expr)
        if(!missing_nrow)
            n <- nrow
        if(missing_ncol)
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
#' @param ... (Optional) Addition \code{axis} argument, specifying the dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{axis = 1}.
#' @return An \linkS4class{Expression} representing the \code{k}th order difference.
#' @examples
#' ## Problem data
#' m <- 101
#' L <- 2
#' h <- L/(m-1)
#'
#' ## Form objective and constraints
#' x <- Variable(m)
#' y <- Variable(m)
#' obj <- sum(y)
#' constr <- list(x[1] == 0, y[1] == 1, x[m] == 1, y[m] == 1, diff(x)^2 + diff(y)^2 <= h^2)
#'
#' ## Solve the catenary problem
#' prob <- Problem(Minimize(obj), constr)
#' result <- solve(prob)
#'
#' ## Plot and compare with ideal catenary
#' xs <- result$getValue(x)
#' ys <- result$getValue(y)
#' plot(c(0, 1), c(0, 1), type = 'n', xlab = "x", ylab = "y")
#' lines(xs, ys, col = "blue", lwd = 2)
#' grid()
#' @docType methods
#' @aliases diff
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
#' @param FUN Hardwired to "*" for the kronecker product.
#' @param make.dimnames (Unimplemented) Dimension names are not supported in \linkS4class{Expression} objects.
#' @param ... (Unimplemented) Optional arguments.
#' @return An \linkS4class{Expression} that represents the kronecker product.
#' @examples
#' X <- cbind(c(1,2), c(3,4))
#' Y <- Variable(2,2)
#' val <- cbind(c(5,6), c(7,8))
#'
#' obj <- X %x% Y
#' prob <- Problem(Minimize(kronecker(X,Y)[1,1]), list(Y == val))
#' result <- solve(prob)
#' result$value
#' result$getValue(kronecker(X,Y))
#' @docType methods
#' @aliases kronecker %x%
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


#'
#' Complex Numbers
#'
#' Basic atoms that support complex arithmetic.
#' 
#' @param z An \linkS4class{Expression} object.
#' @return An \linkS4class{Expression} object that represents the real, imaginary, or complex conjugate.
#' @name complex-atoms
NULL

#' @rdname complex-atoms
#' @export
setMethod("Re", "Expression", function(z) { Real(z) })

#' @rdname complex-atoms
#' @export
setMethod("Im", "Expression", function(z) { Imag(z) })

#' @rdname complex-atoms
#' @export
setMethod("Conj", "Expression", function(z) { if(is_real(z)) z else Conjugate(z) })

#'
#' Product of Entries
#'
#' The product of entries in a vector or matrix.
#'
#' @param expr An \linkS4class{Expression}, vector, or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @return An \linkS4class{Expression} representing the product of the entries of the input.
#' @examples
#' x <- Variable(2)
#' prob <- Problem(Minimize(prod_entries(x)), list(t(x) >= matrix(c(1,2), nrow = 1, ncol = 2)))
#' result <- solve(prob, gp = TRUE)
#' result$value
#' result$getValue(x)
#'
#' C <- Variable(3,2)
#' prob <- Problem(Maximize(prod_entries(C)), list(C[2:3,] <= 2, C[1,] == 1))
#' result <- solve(prob, gp = TRUE)
#' result$value
#' result$getValue(C)
#' @docType methods
#' @name prod_entries
#' @aliases prod
#' @rdname prod_entries
#' @export
prod_entries <- ProdEntries
