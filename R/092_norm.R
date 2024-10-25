## CVXPY SOURCE: cvxpy/atoms/norm.py
#'
#' The Norm atom.
#'
#' Wrapper around the different norm atoms.
#'
#' @param x The matrix to take the norm of
#' @param p The type of norm. Valid options include any positive integer, 'fro' (for frobenius),
#' 'nuc' (sum of singular values), np.inf or 'inf' (infinity norm).
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @return Returns the specified norm of x.
#' @rdname Norm-atom
Norm <- function(x, p = 2, axis = NA_real_, keepdims = FALSE) {
  x <- as.Constant(x)

  # Matrix norms take precedence.
  num_nontrivial_idxs <- sum(dim(x) > 1)
  if(is.na(axis) && ndim(x) == 2) {
    if(p == 1)   # Matrix 1-norm.
      MaxEntries(Norm1(x, axis = 2))
    else if(p == "fro" || (p == 2 && num_nontrivial_idxs == 1))   # Frobenius norm.
      Pnorm(Vec(x), 2)
    else if(p == 2)   # Matrix 2-norm is largest singular value.
      SigmaMax(x)
    else if(p == "nuc")   # The nuclear norm (sum of singular values)
      NormNuc(x)
    else if(p %in% c(Inf, "inf", "Inf"))   # The matrix infinity-norm.
      MaxEntries(Norm1(x, axis = 1))
    else
      stop("Unsupported matrix norm.")
  } else {
    if(p == 1 || is_scalar(x))
      Norm1(x, axis = axis, keepdims = keepdims)
    else if(p %in% c(Inf, "inf", "Inf"))
      NormInf(x, axis = axis, keepdims = keepdims)
    else
      Pnorm(x, p, axis = axis, keepdims = keepdims)
  }
}

## 
#'
#' The Norm2 atom.
#'
#' The 2-norm of an expression.
#'
#' @param x An \linkS4class{Expression} object.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @return Returns the 2-norm of x.
#' @rdname Norm2-atom
Norm2 <- function(x, axis = NA_real_, keepdims = FALSE) {
  Pnorm(x, p = 2, axis = axis, keepdims = keepdims)
}

