## CVXPY SOURCE: cvxpy/reductions/complex2real/canonicalizers/param_canon.py
#'
#' Complex canonicalizer for the parameter matrix atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a parameter matrix atom, where the returned
#' variables are the real component and the imaginary component.
Complex2Real.param_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is_real(expr))
    return(list(expr, NULL))
  else if(is_imag(expr))
    return(list(NULL, Im(expr)))
  else
    return(list(Re(expr), Im(expr)))
}
