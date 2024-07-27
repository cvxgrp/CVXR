## CVXPY SOURCE: cvxpy/reductions/complex2real/canonicalizers/equality_canon.py
#'
#' Complex canonicalizer for the equality constraint
#'
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a equality constraint, where the returned variables are the real component and the imaginary component.
Complex2Real.equality_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is.null(imag_args[[1]]) && is.null(imag_args[[2]]))
    return(list(list(copy(expr, real_args)), NULL))

  # Fill in missing args with zeros.
  for(i in seq_len(imag_args)) {
    if(is.null(imag_args[[i]]))
      imag_args[[i]] <- Constant(matrix(0, nrow = nrow(real_args[[i]]), ncol = ncol(real_args[[i]])))
  }

  imag_cons <- list(EqConstraint(imag_args[[1]], imag_args[[2]], constr_id = real2imag[[as.character(expr@id)]]))

  if(is.null(real_args[[1]]) && is.null(real_args[[2]]))
    return(list(NULL, imag_cons))
  else {
    # Fill in missing args with zeros.
    for(i in seq_len(real_args)) {
      if(is.null(real_args[[i]]))
        real_args[[i]] <- Constant(matrix(0, nrow = nrow(imag_args[[i]]), ncol = ncol(imag_args[[i]])))
    }
    return(list(list(copy(expr, real_args)), imag_cons))
  }
}

#'
#' Complex canonicalizer for the zero constraint
#'
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a zero constraint, where the returned variables are the real component and the imaginary component.
Complex2Real.zero_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is.null(imag_args[[1]]))
    return(list(list(copy(expr, real_args)), NULL))

  imag_cons <- list(ZeroConstraint(imag_args[[1]], constr_id = real2imag[[as.character(expr@id)]]))
  if(is.null(real_args[[1]]))
    return(list(NULL, imag_cons))
  else
    return(list(list(copy(expr, real_args)), imag_cons))
}

