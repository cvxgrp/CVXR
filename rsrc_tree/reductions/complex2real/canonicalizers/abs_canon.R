## CVXPY SOURCE: cvxpy/reductions/complex2real/canonicalizers/abs_canon.py

# Atom canonicalizers.
#'
#' Complex canonicalizer for the absolute value atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of the absolute value atom of a complex expression, where the returned
#' variables are its real and imaginary components parsed out.
Complex2Real.abs_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is.null(real_args[[1]]))   # Imaginary
    output <- abs(imag_args[[1]])
  else if(is.null(imag_args[[1]]))   # Real
    output <- abs(real_args[[1]])
  else {   # Complex
    real <- flatten(real_args[[1]])
    imag <- flatten(imag_args[[1]])
    # norms <- p_norm(hstack(real, imag), p = 2, axis = 1)
    norms <- p_norm(vstack(real, imag), p = 2, axis = 2)
    output <- reshape_expr(norms, dim(real_args[[1]]))
  }
  return(list(output, NULL))
}

