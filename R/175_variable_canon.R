## CVXPY SOURCE: cvxpy/reductions/complex2real/canonicalizers/variable_canon.py
#'
#' Complex canonicalizer for the variable atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a variable atom, where the returned
#' variables are the real component and the NULL imaginary component.
Complex2Real.variable_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is_real(expr))   # Purely real
    return(list(expr, NULL))
  else if(is_imag(expr)) {   # Purely imaginary
    # imag <- Variable(dim(expr), id = real2imag[[as.character()]])
    imag <- new("Variable", dim = dim(expr), id = real2imag[[as.character(expr@id)]])
    return(list(NULL, imag))
  } else if(is_complex(expr) && is_hermitian(expr)) {
    n <- nrow(expr)
    real <- Variable(n, n, var_id = expr@id, symmetric = TRUE)
    if(n > 1) {
      imag_var <- Variable(floor(n*(n-1)/2), var_id = real2imag[[as.character(expr@id)]])
      imag_upper_tri <- UpperTri.vec_to_upper_tri(imag_var, strict = TRUE)
      imag <- SkewSymmetricWrap(imag_upper_tri - t(imag_upper_tri))
    } else
      imag <- Constant(matrix(0))
    return(list(real, imag))
  } else {   # General complex
    expr_dim <- dim(expr)
    real <- new("Variable", dim = expr_dim, var_id = expr@id)
    imag <- new("Variable", dim = expr_dim, var_id = expr@id)
    return(list(real, imag))
  }
}
