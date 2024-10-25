## CVXPY SOURCE: cvxpy/reductions/complex2real/canonicalizers/psd_canon.py
#'
#' Complex canonicalizer for the positive semidefinite atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a positive semidefinite atom, where the returned
#' variables are the real component and the NULL imaginary component.
Complex2Real.psd_canon <- function(expr, real_args, imag_args, real2imag) {
  # Canonicalize functions that take a Hermitian matrix.
  if(is.null(imag_args[[1]]))
    mat <- real_args[[1]]
  else {
    if(is.null(real_args[[1]]))
      real_args[[1]] <- matrix(0, nrow = nrow(imag_args[[1]]), ncol = ncol(imag_args[[1]]))
    mat <- bmat(list(list(real_args[[1]], -imag_args[[1]]),
                     list(imag_args[[1]], real_args[[1]])))
  }
  return(list(list(copy(expr, list(mat))), NULL))
}

