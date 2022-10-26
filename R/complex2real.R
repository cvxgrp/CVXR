
#' 
#' Lifts complex numbers to a real representation.
#' 
#' This reduction takes in a complex problem and returns
#' an equivalent real problem.
#' @rdname Complex2Real-class
Complex2Real <- setClass("Complex2Real", contains = "Reduction")

Complex2Real.accepts <- function(problem) {
  leaves <- c(variables(problem), parameters(problem), constants(problem))
  any(sapply(leaves, function(l) { is_complex(l) }))
}

#' @param object A \linkS4class{Complex2Real} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn Complex2Real Checks whether or not the problem involves any complex numbers.
setMethod("accepts", signature(object = "Complex2Real", problem = "Problem"), function(object, problem) {
  Complex2Real.accepts(problem)
})

#' @describeIn Complex2Real Converts a Complex problem into a Real one.
setMethod("perform", signature(object = "Complex2Real", problem = "Problem"), function(object, problem) {
  inverse_data <- InverseData(problem)

  leaf_map <- list()
  obj <- Complex2Real.canonicalize_tree(problem@objective, inverse_data@real2imag, leaf_map)
  real_obj <- obj[[1]]
  imag_obj <- obj[[2]]

  if(length(imag_obj) > 0)
    stop("Cannot have imaginary component in canonicalized objective")

  constrs <- list()
  for(constraint in problem@constraints) {
    if(inherits(constraint, "EqConstraint"))
      constraint <- lower_equality(constraint)
    constr <- Complex2Real.canonicalize_tree(constraint, inverse_data@real2imag, leaf_map)
    real_constr <- constr[[1]]
    imag_constr <- constr[[2]]
    if(!is.null(real_constr))
      constrs <- c(constrs, real_constr)
    if(!is.null(imag_constr))
      constrs <- c(constrs, imag_constr)
  }

  new_problem <- Problem(real_obj, constrs)
  return(list(object, new_problem, inverse_data))
})

#' @param solution A \linkS4class{Solution} object to invert.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn Complex2Real Returns a solution to the original problem given the inverse data.
setMethod("invert", signature(object = "Complex2Real", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  pvars <- list()
  dvars <- list()
  if(solution@status %in% SOLUTION_PRESENT) {
    for(vid in names(inverse_data@id2var)) {
      var <- inverse_data@id2var[[vid]]
      if(is_real(var))
        pvars[[vid]] <- solution@primal_vars[[vid]]
      else if(is_imag(var)) {
        imag_id <- inverse_data@real2imag[[vid]]
        pvars[[vid]] <- 1i*solution@primal_vars[[as.character(imag_id)]]
      } else if(is_complex(var) && is_hermitian(var)) {
        imag_id <- inverse_data@real2imag[[vid]]
        
        # Imaginary part may have been lost.
        if(as.character(imag_id) %in% names(solution@primal_vars)) {
          imag_val <- solution@primal_vars[[as.character(imag_id)]]
          pvars[[vid]] <- solution@primal_vars[[vid]] + 1i*(imag_val - t(imag_val))/2
        } else
          pvars[[vid]] <- solution@primal_vars[[vid]]
      } else if(is_complex(var)) {
        imag_id <- inverse_data@real2imag[[vid]]
        pvars[[vid]] <- solution@primal_vars[[vid]] + 1i*solution@primal_vars[[as.character(imag_id)]]
      }
    }

    for(cid in names(inverse_data@id2cons)) {
      cons <- inverse_data@id2cons[[cid]]
      if(is_real(cons))
        dvars[[vid]] <- solution@dual_vars[[cid]]
      else if(is_imag(cons)) {
        imag_id <- inverse_data@real2imag[[cid]]
        dvars[[cid]] <- 1i*solution@dual_vars[[as.character(imag_id)]]
      # For equality and inequality constraints.
      } else if((is(cons, "ZeroConstraint") || is(cons, "EqConstraint") || is(cons, "NonPosConstraint")) && is_complex(cons)) {
        imag_id <- inverse_data@real2imag[[cid]]
        dvars[[cid]] <- solution@dual_vars[[cid]] + 1i*solution@dual_vars[[as.character(imag_id)]]
      # For PSD constraints.
      } else if(is(cons, "PSDConstraint") && is_complex(cons)) {
        n <- nrow(cons@args[[1]])
        dual <- solution@dual_vars[[cid]]
        dvars[[cid]] <- dual[1:n,1:n] + 1i*dual[(n+1):nrow(dual), (n+1):ncol(dual)]
      } else
        stop("Unknown constraint type")
    }
  }
  return(Solution(solution@status, solution@opt_val, pvars, dvars, solution@attr))
})

#' 
#' Recursively Canonicalizes a Complex Expression.
#' 
#' @param expr An \linkS4class{Expression} object.
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @param leaf_map A map that consists of a tree representation of the expression.
#' @return A list of the parsed out real and imaginary components of the
#' expression that was constructed by performing the canonicalization of each leaf
#' in the tree.
Complex2Real.canonicalize_tree <- function(expr, real2imag, leaf_map) {
  # TODO: Don't copy affine expressions?
  if(inherits(expr, "PartialProblem"))
    stop("Unimplemented")
  else {
    real_args <- list()
    imag_args <- list()
    for(arg in expr@args) {
      canon <- Complex2Real.canonicalize_tree(arg, real2imag, leaf_map)
      real_args <- c(real_args, list(canon[[1]]))
      imag_args <- c(imag_args, list(canon[[2]]))
    }
    outs <- Complex2Real.canonicalize_expr(expr, real_args, imag_args, real2imag, leaf_map)
    return(list(real_out = outs[[1]], imag_out = outs[[2]]))
  }
}

#' 
#' Canonicalizes a Complex Expression
#' 
#' @param expr An \linkS4class{Expression} object.
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression.
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression.
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @param leaf_map A map that consists of a tree representation of the overall expression
#' @return A list of the parsed out real and imaginary components of the expression at hand.
Complex2Real.canonicalize_expr <- function(expr, real_args, imag_args, real2imag, leaf_map) {
  if(class(expr) %in% names(Complex2Real.CANON_METHODS)) {
    expr_id <- as.character(id(expr))
    # Only canonicalize a variable/constant/parameter once.
    if(length(expr@args) == 0 && expr_id %in% names(leaf_map))
      return(leaf_map[[expr_id]])
    result <- Complex2Real.CANON_METHODS[[class(expr)]](expr, real_args, imag_args, real2imag)
    if(length(expr@args) == 0)
      leaf_map[[expr_id]] <- result
    return(result)
  } else {
    if(!all(sapply(imag_args, is.null)))
      stop("Not all imaginary arguments are NULL")
    return(list(copy(expr, real_args), NULL))
  }
}

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
    norms <- p_norm(hstack(real, imag), p = 2, axis = 1)
    output <- reshape_expr(norms, dim(real_args[[1]]))
  }
  return(list(output, NULL))
}

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
    return(list(0, NULL))
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
    return(list(0, NULL))
  else
    return(list(imag_args[[1]], NULL))
}

#' 
#' Complex canonicalizer for the conjugate atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a conjugate atom, where the returned
#' variables are the real components and negative of the imaginary component.
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
  # 
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
#' @return A canonicalization of a binary atom, where the returned
#' variables are the real component and the imaginary component.
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

#' 
#' Complex canonicalizer for the constant atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a constant atom, where the returned
#' variables are the real component and the imaginary component in the \linkS4class{Constant}
#' atom.
Complex2Real.constant_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is_real(expr))
    return(list(Constant(Re(value(expr))), NULL))
  else if(is_imag(expr))
    return(list(NULL, Constant(Im(value(expr)))))
  else
    return(list(Constant(Re(value(expr))), Constant(Im(value(expr)))))
}

# Matrix canonicalization.
# We expand the matrix A to B = [[Re(A), -Im(A)], [Im(A), Re(A)]]
# B has the same eigenvalues as A (if A is Hermitian).
# If x is an eigenvector of A, then [Re(x), Im(x)] and [Im(x), -Re(x)]
# are eigenvectors with same eigenvalue.
# Thus each eigenvalue is repeated twice.
#' 
#' Complex canonicalizer for the hermitian atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a hermitian matrix atom, where the returned
#' variables are the real component and the imaginary component.
Complex2Real.hermitian_canon <- function(expr, real_args, imag_args, real2imag) {
  # Canonicalize functions that take a Hermitian matrix.
  if(is.null(imag_args[[1]]))
    mat <- real_args[[1]]
  else {
    if(is.null(real_args[[1]]))
      real_args[[1]] <- matrix(0, nrow = nrow(imag_args[[1]]), ncol = ncol(imag_args[[1]]))
    mat <- bmat(list(list(real_args[[1]], -imag_args[[1]]),
                     list(imag_args[[1]], real_args[[1]])))
  }
  return(list(copy(expr, list(mat)), NULL))
}

#' 
#' Complex canonicalizer for the nuclear norm atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a nuclear norm matrix atom, where the returned
#' variables are the real component and the imaginary component.
Complex2Real.norm_nuc_canon <- function(expr, real_args, imag_args, real2imag) {
  # Canonicalize nuclear norm with Hermitian matrix input.
  # Divide by two because each eigenvalue is repeated twice.
  canon <- Complex2Real.hermitian_canon(expr, real_args, imag_args, real2imag)
  real <- canon[[1]]
  imag <- canon[[2]]
  if(!is.null(imag_args[[1]]))
    real <- real/2
  return(list(real, imag))
}

#' 
#' Complex canonicalizer for the largest sum atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of the largest sum atom, where the returned
#' variables are the real component and the imaginary component.
Complex2Real.lambda_sum_largest_canon <- function(expr, real_args, imag_args, real2imag) {
  # Canonicalize nuclear norm with Hermitian matrix input.
  # Divide by two because each eigenvalue is repeated twice.
  canon <- Complex2Real.hermitian_canon(expr, real_args, imag_args, real2imag)
  real <- canon[[1]]
  imag <- canon[[2]]
  real@k <- 2*real@k
  if(!is.null(imag_args[[1]]))
    real <- real/2
  return(list(real, imag))
}

#'
#' Upcast 0D and 1D to 2D.
#'
#' @param expr An \linkS4class{Expression} object
#' @return An expression of dimension at least 2. 
Complex2Real.at_least_2D <- function(expr) {
  # 
  if(length(dim(expr)) < 2)
    return(reshape_expr(expr, c(size(expr), 1)))
  else
    return(expr)
}

#' 
#' Complex canonicalizer for the quadratic atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a quadratic atom, where the returned
#' variables are the real component and the imaginary component as NULL.
Complex2Real.quad_canon <- function(expr, real_args, imag_args, real2imag) {
  # Convert quad_form to real.
  if(is.null(imag_args[[1]])) {
    vec <- real_args[[1]]
    mat <- real_args[[2]]
  } else if(is.null(real_args[[1]])) {
    vec <- imag_args[[1]]
    mat <- real_args[[2]]
  } else {
    vec <- vstack(Complex2Real.at_least_2D(real_args[[1]]),
                  Complex2Real.at_least_2D(imag_args[[1]]))
    if(is.null(real_args[[2]]))
      real_args[[2]] <- matrix(0, nrow = nrow(imag_args[[2]]), ncol = ncol(imag_args[[2]]))
    else if(is.null(imag_args[[2]]))
      imag_args[[2]] <- matrix(0, nrow = nrow(real_args[[2]]), ncol = ncol(real_args[[2]]))
    mat <- bmat(list(list(real_args[[2]], -imag_args[[2]]),
                     list(imag_args[[2]], real_args[[2]])
                ))
    # mat <- Constant(value(mat))
    mat <- PSDWrap(mat)
  }
  return(list(copy(expr, list(vec, mat)), NULL))
}

#' 
#' Complex canonicalizer for the quadratic over linear term atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a quadratic over a linear term atom, where the returned
#' variables are the real component and the imaginary component.
Complex2Real.quad_over_lin_canon <- function(expr, real_args, imag_args, real2imag) {
  # Convert quad_over_lin to real.
  mat <- bmat(list(real_args[[1]], imag_args[[1]]))
  return(list(copy(expr, list(mat, real_args[[2]])), NULL))
}


#'
#' Complex canonicalizer for the matrix fraction atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a matrix atom, where the returned
#' variables are converted to real variables.
Complex2Real.matrix_frac_canon <- function(expr, real_args, imag_args, real2imag) {
  # Convert matrix_frac to real.
  if(is.null(real_args[[1]]))
    real_args[[1]] <- matrix(0, nrow = nrow(imag_args[[1]]), ncol = ncol(imag_args[[1]]))
  if(is.null(imag_args[[1]]))
    imag_args[[1]] <- matrix(0, nrow = nrow(real_args[[1]]), ncol = ncol(real_args[[1]]))
  vec <- vstack(Complex2Real.at_least_2D(real_args[[1]]),
                Complex2Real.at_least_2D(imag_args[[1]]))
  if(is.null(real_args[[2]]))
    real_args[[2]] <- matrix(0, nrow = nrow(imag_args[[2]]), ncol = ncol(imag_args[[2]]))
  else if(is.null(imag_args[[2]]))
    imag_args[[2]] <- matrix(0, nrow = nrow(real_args[[2]]), ncol = ncol(real_args[[2]]))
  mat <- bmat(list(list(real_args[[2]], -imag_args[[2]]),
                   list(imag_args[[2]], real_args[[2]])
              ))
  return(list(copy(expr, list(vec, mat)), NULL))
}

#' 
#' Complex canonicalizer for the non-positive atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a non positive atom, where the returned
#' variables are the real component and the imaginary component.
Complex2Real.nonpos_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is.null(imag_args[[1]]))
    return(list(list(copy(expr, real_args)), NULL))
  
  imag_cons <- list(NonPosConstraint(imag_args[[1]], id = real2imag[[as.character(id(expr))]]))
  if(is.null(real_args[[1]]))
    return(list(NULL, imag_cons))
  else
    return(list(list(copy(expr, real_args)), imag_cons))
}

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
  else if(is_imag(expr)) {
    imag <- CallbackParam(function() { list(Im(value(expr)), dim(expr)) })
    return(list(NULL, imag))
  } else {
    real <- CallbackParam(function() { list(Re(value(expr)), dim(expr)) })
    imag <- CallbackParam(function() { list(Im(value(expr)), dim(expr)) })
    return(list(real, imag))
  }
}

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

#' 
#' Complex canonicalizer for the SOC atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a SOC atom, where the returned
#' variables are the real component and the NULL imaginary component.
Complex2Real.soc_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is.null(real_args[[2]]))   # Imaginary.
    output <- list(SOC(real_args[[1]], imag_args[[2]], axis = expr@axis, id = real2imag[[as.character(expr@id)]]))
  else if(is.null(imag_args[[2]]))   # Real.
    output <- list(SOC(real_args[[1]], real_args[[2]], axis = expr@axis, id = expr@id))
  else {   # Complex.
    orig_dim <- dim(real_args[[2]])
    real <- flatten(real_args[[2]])
    imag <- flatten(imag_args[[2]])
    flat_X <- new("Variable", dim = dim(real))
    inner_SOC <- SOC(flat_X, vstack(real, imag), axis = 1)   # TODO: Check the axis here is correct.
    real_X <- reshape_expr(flat_X, orig_dim)
    outer_SOC <- SOC(real_args[[1]], real_X, axis = expr@axis, id = expr@id)
    output <- list(inner_SOC, outer_SOC)
  }
  return(list(output, NULL))
}

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
  if(is_real(expr))
    return(list(expr, NULL))

  # imag <- Variable(dim(expr), id = real2imag[[as.character(id(expr))]])
  imag <- new("Variable", dim = dim(expr), id = real2imag[[as.character(expr@id)]])
  if(is_imag(expr))
    return(list(NULL, imag))
  else if(is_complex(expr) && is_hermitian(expr))
    # return(list(Variable(dim(expr), id = id(expr), symmetric = TRUE), (imag - t(imag))/2))
    return(list(new("Variable", dim = dim(expr), id = expr@id, symmetric = TRUE), (imag - t(imag))/2))
  else   # Complex.
    # return(list(Variable(dim(expr), id = id(expr)), imag))
    return(list(new("Variable", dim = dim(expr), id = expr@id), imag))
}

#' 
#' Complex canonicalizer for the zero atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param real_args A list of \linkS4class{Constraint} objects for the real part of the expression
#' @param imag_args A list of \linkS4class{Constraint} objects for the imaginary part of the expression
#' @param real2imag A list mapping the ID of the real part of a complex expression to the ID of its imaginary part.
#' @return A canonicalization of a zero atom, where the returned
#' variables are the real component and the imaginary component.
Complex2Real.zero_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is.null(imag_args[[1]]))
    return(list(list(copy(expr, real_args)), NULL))

  imag_cons <- list(ZeroConstraint(imag_args[[1]], id = real2imag[[as.character(id(expr))]]))
  if(is.null(real_args[[1]]))
    return(list(NULL, imag_cons))
  else
    return(list(list(copy(expr, real_args)), imag_cons))
}

Complex2Real.CANON_METHODS <- list(AddExpression = Complex2Real.separable_canon,
                                   Bmat = Complex2Real.separable_canon,
                                   CumSum = Complex2Real.separable_canon,
                                   Diag = Complex2Real.separable_canon,
                                   HStack = Complex2Real.separable_canon,
                                   Index = Complex2Real.separable_canon,
                                   SpecialIndex = Complex2Real.separable_canon,
                                   Promote = Complex2Real.separable_canon,
                                   Reshape = Complex2Real.separable_canon,
                                   SumEntries = Complex2Real.separable_canon,
                                   Trace = Complex2Real.separable_canon,
                                   Transpose = Complex2Real.separable_canon,
                                   NegExpression = Complex2Real.separable_canon,
                                   UpperTri = Complex2Real.separable_canon,
                                   VStack = Complex2Real.separable_canon,
                                   
                                   Conv = Complex2Real.binary_canon,
                                   DivExpression = Complex2Real.binary_canon,
                                   Kron = Complex2Real.binary_canon,
                                   MulExpression = Complex2Real.binary_canon,
                                   Multiply = Complex2Real.binary_canon,
                                   
                                   Conjugate = Complex2Real.conj_canon,
                                   Imag = Complex2Real.imag_canon,
                                   Real = Complex2Real.real_canon,
                                   Variable = Complex2Real.variable_canon,
                                   Constant = Complex2Real.constant_canon,
                                   Parameter = Complex2Real.param_canon,
                                   NonPosConstraint = Complex2Real.nonpos_canon,
                                   PSDConstraint = Complex2Real.psd_canon,
                                   SOC = Complex2Real.soc_canon,
                                   ZeroConstraint = Complex2Real.zero_canon,
                                   
                                   Abs = Complex2Real.abs_canon,
                                   Norm1 = Complex2Real.pnorm_canon,
                                   NormInf = Complex2Real.pnorm_canon,
                                   Pnorm = Complex2Real.pnorm_canon,
                                   
                                   LambdaMax = Complex2Real.hermitian_canon,
                                   LogDet = Complex2Real.norm_nuc_canon,
                                   NormNuc = Complex2Real.norm_nuc_canon,
                                   SigmaMax = Complex2Real.hermitian_canon,
                                   QuadForm = Complex2Real.quad_canon,
                                   QuadOverLin = Complex2Real.quad_over_lin_canon,
                                   MatrixFrac = Complex2Real.matrix_frac_canon,
                                   LambdaSumLargest = Complex2Real.lambda_sum_largest_canon)
