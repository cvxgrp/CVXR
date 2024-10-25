## CVXPY SOURCE: cvxpy/atoms/quad_form.py

#'
#' The QuadForm class.
#'
#' This class represents the quadratic form \eqn{x^T P x}
#'
#' @slot x An \linkS4class{Expression} or numeric vector.
#' @slot P An \linkS4class{Expression}, numeric matrix, or vector.
#' @name QuadForm-class
#' @aliases QuadForm
#' @rdname QuadForm-class
.QuadForm <- setClass("QuadForm", representation(x = "ConstValORExpr", P = "ConstValORExpr"), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric vector.
#' @param P An \linkS4class{Expression}, numeric matrix, or vector.
#' @rdname QuadForm-class
QuadForm <- function(x, P) { .QuadForm(x = x, P = P) }

setMethod("initialize", "QuadForm", function(.Object, ..., x, P) {
  .Object@x <- x
  .Object@P <- P
  callNextMethod(.Object, ..., atom_args = list(.Object@x, .Object@P))
})

#' @describeIn QuadForm The name and arguments of the atom.
setMethod("name", "QuadForm", function(x) {
  paste(class(x), "(", x@args[[1]], ", ", x@args[[2]], ")", sep = "")
})

#' @param object A \linkS4class{QuadForm} object.
#' @describeIn QuadForm Does the atom handle complex numbers?
setMethod("allow_complex", "QuadForm", function(object) { TRUE })

#' @param values A list of numeric values for the arguments
#' @describeIn QuadForm Returns the quadratic form.
setMethod("to_numeric", "QuadForm", function(object, values) {
  prod <- values[[2]] %*% values[[1]]
  if(is_complex(object@args[[1]]))
    return(t(Conj(values[[1]])) %*% prod)
  else
    return(t(values[[1]]) %*% prod)
})

#' @describeIn QuadForm Checks the dimensions of the arguments.
setMethod("validate_args", "QuadForm", function(object) {
  callNextMethod()
  n <- nrow(object@args[[2]])
  x_dim <- dim(object@args[[1]])
  # if(ncol(object@args[[2]]) != n || !(dim(object@args[[1]]) %in% list(c(n, 1), c(n, NA_real_))))
  if(ncol(object@args[[2]]) != n || !(length(x_dim) == 2 && all(x_dim == c(n,1))))
    stop("Invalid dimensions for arguments.")
})

#' @describeIn QuadForm Returns the sign (is positive, is negative) of the atom.
setMethod("sign_from_args", "QuadForm", function(object) { c(is_atom_convex(object), is_atom_concave(object)) })

#' @describeIn QuadForm The dimensions of the atom.
setMethod("dim_from_args", "QuadForm", function(object) {
  # if(ndim(object@args[[1]]) == 0)
  #  c()
  # else
  #  c(1,1)
  c(1,1)
})

#' @describeIn QuadForm Is the atom convex?
setMethod("is_atom_convex", "QuadForm", function(object) { is_psd(object@args[[2]]) })

#' @describeIn QuadForm Is the atom concave?
setMethod("is_atom_concave", "QuadForm", function(object) { is_nsd(object@args[[2]]) })

#' @describeIn QuadForm Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "QuadForm", function(object) { TRUE })

#' @describeIn QuadForm Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "QuadForm", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn QuadForm Is the atom weakly increasing in the argument \code{idx}?
setMethod("is_incr", "QuadForm", function(object, idx) {
  (is_nonneg(object@args[[1]]) && is_nonneg(object@args[[2]])) ||
    (is_nonpos(object@args[[1]]) && is_nonneg(object@args[[2]]))
})

#' @param idx An index into the atom.
#' @describeIn QuadForm Is the atom weakly decreasing in the argument \code{idx}?
setMethod("is_decr", "QuadForm", function(object, idx) {
  (is_nonneg(object@args[[1]]) && is_nonpos(object@args[[2]])) ||
    (is_nonpos(object@args[[1]]) && is_nonpos(object@args[[2]]))
})

#' @describeIn QuadForm Is the atom quadratic?
setMethod("is_quadratic", "QuadForm", function(object) { TRUE })

#' @describeIn QuadForm Is the atom piecewise linear?
setMethod("is_pwl", "QuadForm", function(object) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn QuadForm Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "QuadForm", function(object, values) {
  x <- values[[1]]
  P <- values[[2]]
  D <- 2*P %*% t(x)
  Matrix(as.vector(t(D)), sparse = TRUE)
})

#'
#' Compute a Matrix Decomposition.
#'
#' Compute sgn, scale, M such that \eqn{P = sgn * scale * dot(M, t(M))}.
#'
#' @param P A real symmetric positive or negative (semi)definite input matrix
#' @param cond Cutoff for small eigenvalues. Singular values smaller than rcond * largest_eigenvalue are considered negligible.
#' @param rcond Cutoff for small eigenvalues. Singular values smaller than rcond * largest_eigenvalue are considered negligible.
#' @return A list consisting of induced matrix 2-norm of P and a rectangular matrix such that P = scale * (dot(M1, t(M1)) - dot(M2, t(M2)))
.decomp_quad <- function(P, cond = NA, rcond = NA) {
  eig <- eigen(P, only.values = FALSE)
  w <- eig$values
  V <- eig$vectors

  if(!is.na(rcond))
    cond <- rcond
  if(cond == -1 || is.na(cond))
    cond <- 1e6 * .Machine$double.eps   # All real numbers are stored as double precision in R

  scale <- max(abs(w))
  if(scale < cond)
    return(list(scale = 0, M1 = V[,FALSE], M2 = V[,FALSE]))
  w_scaled <- w / scale
  maskp <- w_scaled > cond
  maskn <- w_scaled < -cond

  # TODO: Allow indefinite QuadForm
  if(any(maskp) && any(maskn))
    warning("Forming a non-convex expression QuadForm(x, indefinite)")

  if(sum(maskp) <= 1)
    M1 <- as.matrix(V[,maskp] * sqrt(w_scaled[maskp]))
  else
    M1 <- V[,maskp] %*% diag(sqrt(w_scaled[maskp]))

  if(sum(maskn) <= 1)
    M2 <- as.matrix(V[,maskn]) * sqrt(-w_scaled[maskn])
  else
    M2 <- V[,maskn] %*% diag(sqrt(-w_scaled[maskn]))
  list(scale = scale, M1 = M1, M2 = M2)
}

#'
#' The SymbolicQuadForm class.
#'
#' @slot x An \linkS4class{Expression} or numeric vector.
#' @slot P An \linkS4class{Expression}, numeric matrix, or vector.
#' @slot original_expression The original \linkS4class{Expression}.
#' @name SymbolicQuadForm-class
#' @aliases SymbolicQuadForm
#' @rdname SymbolicQuadForm-class
.SymbolicQuadForm <- setClass("SymbolicQuadForm", representation(x = "ConstValORExpr", P = "ConstValORExpr", original_expression = "Expression"), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric vector.
#' @param P An \linkS4class{Expression}, numeric matrix, or vector.
#' @param expr The original \linkS4class{Expression}.
#' @rdname SymbolicQuadForm-class
SymbolicQuadForm <- function(x, P, expr) { .SymbolicQuadForm(x = x, P = P, original_expression = expr) }

setMethod("initialize", "SymbolicQuadForm", function(.Object, ..., x, P, original_expression) {
  .Object@x <- x
  .Object@original_expression <- original_expression
  .Object <- callNextMethod(.Object, ..., atom_args = list(x, P), validate = FALSE)
  .Object@P <- .Object@args[[2]]
  validObject(.Object)
  .Object
})

#' @param object A \linkS4class{SymbolicQuadForm} object.
#' @describeIn SymbolicQuadForm The dimensions of the atom.
setMethod("dim_from_args", "SymbolicQuadForm", function(object) { dim_from_args(object@original_expression) })

#' @describeIn SymbolicQuadForm The sign (is positive, is negative) of the atom.
setMethod("sign_from_args", "SymbolicQuadForm", function(object) { sign_from_args(object@original_expression) })

#' @describeIn SymbolicQuadForm The original expression.
setMethod("get_data", "SymbolicQuadForm", function(object) { list(object@original_expression) })

#' @describeIn SymbolicQuadForm Is the original expression convex?
setMethod("is_atom_convex", "SymbolicQuadForm", function(object) { is_atom_convex(object@original_expression) })

#' @describeIn SymbolicQuadForm Is the original expression concave?
setMethod("is_atom_concave", "SymbolicQuadForm", function(object) { is_atom_concave(object@original_expression) })

#' @param idx An index into the atom.
#' @describeIn SymbolicQuadForm Is the original expression weakly increasing in argument \code{idx}?
setMethod("is_incr", "SymbolicQuadForm", function(object, idx) { is_incr(object@original_expression, idx) })

#' @describeIn SymbolicQuadForm Is the original expression weakly decreasing in argument \code{idx}?
setMethod("is_decr", "SymbolicQuadForm", function(object, idx) { is_decr(object@original_expression, idx) })

#' @describeIn SymbolicQuadForm The atom is quadratic.
setMethod("is_quadratic", "SymbolicQuadForm", function(object) { TRUE })

#' @param values A list of numeric values for the arguments
#' @describeIn SymbolicQuadForm Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "SymbolicQuadForm", function(object, values) { stop("Unimplemented") })

