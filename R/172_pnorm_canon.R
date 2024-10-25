## CVXPY SOURCE: cvxpy/reductions/complex2real/canonicalizers/pnorm_canon.py

#'
#' Complex canonicalizer for the p norm atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a pnorm atom, where the returned
#' variables are the real component and the NULL imaginary component.
Complex2Real.pnorm_canon <- function(expr, real_args, imag_args, real2imag) {
  abs_args <- Complex2Real.abs_canon(expr, real_args, imag_args, real2imag)
  abs_real_args <- abs_args[[1]]
  return(list(copy(expr, list(abs_real_args)), NULL))
}

