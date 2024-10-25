## CVXPY SOURCE: cvxpy/reductions/complex2real/canonicalizers/aff_canon.py
# Affine canonicalization.
#'
#' Complex canonicalizer for the separable atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a separable atom, where the returned
#' variables are its real and imaginary components parsed out.
Complex2Real.separable_canon <- function(expr, real_args, imag_args, real2imag) {
  # Canonicalize linear functions that are separable in real and imaginary parts.
  if(all(sapply(imag_args, is.null)))
    outputs <- list(copy(expr, real_args), NULL)
  else if(all(sapply(real_args, is.null)))
    outputs <- list(NULL, copy(expr, imag_args))
  else {   # Mixed real and imaginary arguments.
    for(idx in seq_along(real_args)) {
      real_val <- real_args[[idx]]
      if(is.null(real_val))
        real_args[[idx]] <- Constant(matrix(0, nrow = nrow(imag_args[[idx]]), ncol = ncol(imag_args[[idx]])))
      else if(is.null(imag_args[[idx]]))
        imag_args[[idx]] <- Constant(matrix(0, nrow = nrow(real_args[[idx]]), ncol = ncol(real_args[[idx]])))
    }
    outputs <- list(copy(expr, real_args), copy(expr, imag_args))
  }
  return(outputs)
}

#'
#' Complex canonicalizer for the real atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a real atom, where the returned
#' variables are the real component and NULL for the imaginary component.
Complex2Real.real_canon <- function(expr, real_args, imag_args, real2imag) {
  # If no real arguments, return zero.
  if(is.null(real_args[[1]]))
    return(list(0*imag_args[[1]], NULL))
  else
    return(list(real_args[[1]], NULL))
}

#'
#' Complex canonicalizer for the imaginary atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of an imaginary atom, where the returned
#' variables are the imaginary component and NULL for the real component.
Complex2Real.imag_canon <- function(expr, real_args, imag_args, real2imag) {
  # If no imaginary arguments, return zero.
  if(is.null(imag_args[[1]]))
    return(list(0*real_args[[1]], NULL))
  else
    return(list(imag_args[[1]], NULL))
}

#'
#' Wrapper for Hermitian canonicalizer
#'
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a hermitian wrapper, where the returned variables are the real components and negative of the imaginary component.
Complex2Real.hermitian_wrap_canon <- function(expr, real_args, imag_args, real2imag) {
  if(!is.null(imag_args[[1]]))
    imag_arg <- SkewSymmetricWrap(imag_args[[1]])
  else {
    # This is a weird code path to hit.
    imag_arg <- NULL
  }
  real_arg <- SymmetricWrap(real_args[[1]])
  return(list(real_arg, imag_arg))
}

#'
#' Complex canonicalizer for the conjugate atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a conjugate atom, where the returned variables are the real components and negative of the imaginary component.
Complex2Real.conj_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is.null(imag_args[[1]]))
    imag_arg <- NULL
  else
    imag_arg <- -imag_args[[1]]
  return(list(real_args[[1]], imag_arg))
}

#'
#' Helper function to combine arguments.
#'
#' @param expr An \linkS4class{Expression} object
#' @param lh_arg The arguments for the left-hand side
#' @param rh_arg The arguments for the right-hand side
#' @return A joined expression of both left and right expressions
Complex2Real.join <- function(expr, lh_arg, rh_arg) {
  if(is.null(lh_arg) || is.null(rh_arg))
    return(NULL)
  else
    return(copy(expr, list(lh_arg, rh_arg)))
}

#'
#' Helper function to sum arguments.
#'
#' @param lh_arg The arguments for the left-hand side
#' @param rh_arg The arguments for the right-hand side
#' @param neg Whether to negate the right hand side
Complex2Real.add <- function(lh_arg, rh_arg, neg = FALSE) {
  # Negates rh_arg if neg is TRUE.
  if(!is.null(rh_arg) && neg)
    rh_arg <- -rh_arg

  if(is.null(lh_arg) && is.null(rh_arg))
    return(NULL)
  else if(is.null(lh_arg))
    return(rh_arg)
  else if(is.null(rh_arg))
    return(lh_arg)
  else
    return(lh_arg + rh_arg)
}

#'
#' Complex canonicalizer for the binary atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a binary atom, where the returned variables are the real component and the imaginary component.
Complex2Real.binary_canon <- function(expr, real_args, imag_args, real2imag) {
  # Canonicalize functions like multiplication.
  real_by_real <- Complex2Real.join(expr, real_args[[1]], real_args[[2]])
  imag_by_imag <- Complex2Real.join(expr, imag_args[[1]], imag_args[[2]])
  real_by_imag <- Complex2Real.join(expr, real_args[[1]], imag_args[[2]])
  imag_by_real <- Complex2Real.join(expr, imag_args[[1]], real_args[[2]])
  real_output <- Complex2Real.add(real_by_real, imag_by_imag, neg = TRUE)
  imag_output <- Complex2Real.add(real_by_imag, imag_by_real, neg = FALSE)
  return(list(real_output, imag_output))
}
