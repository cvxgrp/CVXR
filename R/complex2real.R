# Lifts complex numbers to a real representation.
setClass("Complex2Real", contains = "Reduction")

Complex2Real.accepts <- function(problem) {
  leaves <- c(variables(problem), parameters(problem), constants(problem))
  any(sapply(leaves, function(l) { is_complex(l) }))
}

setMethod("accepts", signature(object = "Complex2Real", problem = "Problem"), function(object, problem) {
  Complex2Real.accepts(problem)
})

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
    if(class(constraint) == "EqConstraint")
      constraint <- lower_equality(constraint)
    constr <- Complex2Real.canonicalize_tree(constraint, inverse_data@real2imag, leaf_map)
    real_constr <- constr[[1]]
    imag_constr <- constr[[2]]
    if(!is.na(real_constr))
      constrs <- c(constrs, real_constr)
    if(!is.na(imag_constr))
      constrs <- c(constrs, imag_constr)
  }

  new_problem <- Problem(real_obj, constrs)
  return(list(new_problem, inverse_data))
})

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
        pvars[[vid]] <- 1i*solution@primal_vars[[vid]]
      } else if(is_complex(var) && is_hermitian(var)) {
        imag_id <- inverse_data@real2imag[[vid]]
        imag_val <- solution@primal_vars[[imag_id]]
        pvars[[vid]] <- solution@primal_vars[[vid]] + 1i*(imag_val - t(imag_val))/2
      } else if(is_complex(var)) {
        imag_id <- inverse_data@real2imag[[vid]]
        pvars[[vid]] <- solution@primal_vars[[vid]] + 1i*solution@primal_vars[[imag_id]]
      }
    }

    for(cid in names(inverse_data@id2cons)) {
      cons <- inverse_data@id2cons[[cid]]
      if(is_real(cons))
        dvars[[vid]] <- solution@dual_vars[[cid]]
      else if(is_imag(cons)) {
        imag_id <- inverse_data@real2imag[[cid]]
        dvars[[cid]] <- 1i*solution@dual_vars[[imag_id]]
      # For equality and inequality constraints.
      } else if((is(cons, "ZeroConstraint") || is(cons, "EqConstraint") || is(cons, "NonPosConstraint")) && is_complex(cons)) {
        imag_id <- inverse_data@real2imag[[cid]]
        dvars[[cid]] <- solution@dual_vars[[cid]] + 1i*solution@dual_vars[[imag_id]]
      # For PSD constraints.
      } else if(is(cons, "PSD") && is_complex(cons)) {
        n <- cons@args[[1]]@dim[1]
        dual <- solution@dual_vars[[cid]]
        dvars[[cid]] <- dual[1:n,1:n] + 1i*dual[(n+1):nrow(dual), (n+1):ncol(dual)]
      } else
        stop("Unknown constraint type")
    }
  }
  return(Solution(solution@status, solution@opt_val, pvars, dvars, solution@attr))
})

Complex2Real.canonicalize_tree <- function(expr, real2imag, leaf_map) {
  # TODO: Don't copy affine expressions?
  if(type(expr) == "PartialProblem")
    stop("Unimplemented")
  else {
    real_args <- list()
    imag_args <- list()
    for(arg in expr@args) {
      canon <- Complex2Real.canonicalize_tree(arg, real2imag, leaf_map)
      real_args <- c(real_args, canon[[1]])
      imag_args <- c(imag_args, canon[[2]])
    }
    outs <- Complex2Real.canonicalize_expr(expr, real_args, imag_args, real2imag, leaf_map)
    return(list(real_out = outs[[1]], imag_out = outs[[2]]))
  }
}

Complex2Real.canonicalize_expr <- function(expr, real_args, imag_args, real2imag, leaf_map) {
  if(class(expr) %in% names(Complex2Real.CANON_METHODS)) {
    # Only canonicalize a variable/constant/parameter once.
    if(length(expr@args) == 0 && expr %in% leaf_map)
      return(leaf_map[expr])
    result <- elim_cplx_methods[type(expr)](expr, real_args, imag_args, real2imag)
    if(length(expr@args) == 0)
      leaf_map[expr] <- result
    return(result)
  } else {
    if(!all(sapply(imag_args, function(v) { is.na(v) || is.null(v) })))
      stop("Not all imaginary arguments are NA or NULL")
    return(list(copy(expr, real_args), NA))
  }
}

# Atom canonicalizers.
Complex2Real.abs_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is.na(real_args[[1]]))   # Imaginary
    output <- abs(imag_args[[1]])
  else if(is.na(imag_args[[1]]))   # Real
    output <- abs(real_args[[1]])
  else {   # Complex
    real <- real_args[[1]]
    imag <- imag_args[[1]]
    norms <- p_norm(vstack(c(real, imag)), p = 2, axis = 2)
    output <- reshape(norms, real_args[[1]]@dim)
  }
  return(list(output, NA))
}

# Affine canonicalization.
Complex2Real.separable_canon <- function(expr, real_args, imag_args, real2imag) {
  # Canonicalize linear functions that are separable in real and imaginary parts.
  if(all(is.na(imag_args)))
    outputs <- list(copy(expr, real_args), NA)
  else if(all(is.na(real_args)))
    ouputs <- list(NA, copy(expr, imag_args))
  else {   # Mixed real and imaginary arguments.
    for(idx in length(real_args)) {
      real_val <- real_args[idx]
      if(is.na(real_val))
        real_args[idx] <- Constant(matrix(0, nrow = nrow(imag_args[idx]), ncol = ncol(imag_args[idx])))
      else if(is.na(imag_args[idx]))
        imag_args[idx] <- Constant(matrix(0, nrow = nrow(real_args[idx]), ncol = ncol(real_args[idx])))
    }
    outputs <- list(copy(expr, real_args), copy(expr, imag_args))
  }
  return(outputs)
}

Complex2Real.real_canon <- function(expr, real_args, imag_args, real2imag) {
  list(real_args[[1]], NA)
}

Complex2Real.imag_canon <- function(expr, real_args, imag_args, real2imag) {
  list(imag_args[[1]], NA)
}

Complex2Real.conj_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is.na(imag_args[[1]]))
    imag_arg <- NA
  else
    imag_arg <- -imag_args[[1]]
  return(list(real_args[[1]], imag_arg))
}

Complex2Real.join <- function(expr, lh_arg, rh_arg) {
  # Helper function to combine arguments.
  if(is.na(lh_arg) || is.na(rh_arg))
    return(NA)
  else
    return(copy(expr, list(lh_arg, rh_arg)))
}

Complex2Real.add <- function(lh_arg, rh_arg, neg = FALSE) {
  # Helper function to sum arguments.
  # Negates rh_arg if neg is TRUE.
  if(!is.na(rh_arg) && neg)
    rh_arg <- -rh_arg

  if(is.na(lh_arg) && is.na(rh_arg))
    return(NA)
  else if(is.na(lh_arg))
    return(rh_arg)
  else if(is.na(rh_arg))
    return(lh_arg)
  else
    return(lh_arg + rh_arg)
}

Complex2Real.binary_canon <- function(expr, real_args, imag_args, real2imag) {
  # Canonicalize functions like multiplication.
  real_by_real <- Complex2Real.join(expr, real_args[[1]], real_args[[2]])
  imag_by_imag <- Complex2Real.join(expr, imag_args[[1]], imag_args[[2]])
  real_by_imag <- Complex2Real.join(expr, real_args[[1]], imag_args[[2]])
  imag_by_real <- Complex2Real.join(expr, imag_args[[1]], real_args[[2]])
  real_output <- Complex2Real.add(real_by_real, imag_by_imag, neg = TRUE)
  imag_output <- Complex2Real.add(real_by_imag, imag_by_real, neg = TRUE)
  return(list(real_output, imag_output))
}

Complex2Real.constant_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is_real(expr))
    return(list(Constant(Re(value(expr))), NA))
  else if(is_imag(expr))
    return(list(NA, Constant(Im(value(expr)))))
  else
    return(list(Constant(Re(value(expr))), Constant(Im(value(expr)))))
}

# Matrix canonicalization.
# We expand the matrix A to B = [[Re(A), -Im(A)], [Im(A), Re(A)]]
# B has the same eigenvalues as A (if A is Hermitian).
# If x is an eigenvector of A, then [Re(x), Im(x)] and [Im(x), -Re(x)]
# are eigenvectors with same eigenvalue.
# Thus each eigenvalue is repeated twice.
Complex2Real.hermitian_canon <- function(expr, real_args, imag_args, real2imag) {
  # Canonicalize functions taht take a Hermitian matrix.
  if(is.na(imag_args[[1]]))
    mat <- real_args[[1]]
  else {
    if(is.na(real_args[[1]]))
      real_args[[1]] <- matrix(0, nrow = nrow(imag_args[[1]]), ncol = ncol(imag_args[[1]]))
    mat <- bmat(list(list(real_args[[1]], -imag_args[[1]]),
                     list(imag_args[[1]], real_args[[1]])
                ))
  }
  return(list(copy(expr, list(mat)), NA))
}

Complex2Real.norm_nuc_canon <- function(expr, real_args, imag_args, real2imag) {
  # Canonicalize nuclear norm with Hermitian matrix input.
  # Divide by two because each eigenvalue is repeated twice.
  canon <- Complex2Real.hermitian_canon(expr, real_args, imag_args, real2imag)
  real <- canon[[1]]
  imag <- canon[[2]]
  if(!is.na(imag_args[[1]]))
    real <- real/2
  return(list(real, imag))
}

Complex2Real.lambda_sum_largest_canon <- function(expr, real_args, imag_args, real2imag) {
  # Canonicalize nuclear norm with Hermitian matrix input.
  # Divide by two because each eigenvalue is repeated twice.
  canon <- Complex2Real.hermitian_canon(expr, real_args, imag_args, real2imag)
  real <- canon[[1]]
  imag <- canon[[2]]
  real@k <- 2*real@k
  if(!is.na(imag_args[[1]]))
    real <- real/2
  return(list(real, imag))
}

Complex2Real.at_least_2D <- function(expr) {
  # Upcast 0D and 1D to 2D.
  if(ndim(expr) < 2)
    return(reshape(expr, c(size(expr), 1)))
  else
    return(expr)
}

Complex2Real.quad_canon <- function(expr, real_args, imag_args, real2imag) {
  # Convert quad_form to real.
  if(is.na(imag_args[[1]])) {
    vec <- real_args[[1]]
    mat <- real_args[[2]]
  } else if(is.na(real_args[[1]])) {
    vec <- imag_args[[1]]
    mat <- real_args[[2]]
  } else {
    vec <- vstack(list(Complex2Real.at_least_2D(real_args[[1]]),
                       Complex2Real.at_least_2D(imag_args[[1]])))
    if(is.na(real_args[[2]]))
      real_args[[2]] <- matrix(0, nrow = nrow(imag_args[[2]]), ncol = ncol(imag_args[[2]]))
    else if(is.na(imag_args[[2]]))
      imag_args[[2]] <- matrix(0, nrow = nrow(real_args[[2]]), ncol = ncol(real_args[[2]]))
    mat <- bmat(list(list(real_args[[2]], -imag_args[[2]]),
                     list(imag_args[[2]], real_args[[2]])
                ))
    # HACK TODO
    mat <- Constant(value(mat))
  }
  return(list(copy(expr, list(vec, mat)), NA))
}

Complex2Real.matrix_frac_canon <- function(expr, real_args, imag_args, real2imag) {
  # Convert matrix_frac to real.
  if(is.na(real_args[[1]]))
    real_args[[1]] <- matrix(0, nrow = nrow(imag_args[[1]]), ncol = ncol(imag_args[[1]]))
  if(is.na(imag_args[[1]]))
    imag_args[[1]] <- matrix(0, nrow = nrow(real_args[[1]]), ncol = ncol(real_args[[1]]))
  vec <- vstack(list(Complex2Real.at_least_2D(real_args[[1]]),
                     Complex2Real.at_least_2D(imag_args[[1]])))
  if(is.na(real_args[[2]]))
    real_args[[2]] <- matrix(0, nrow = nrow(imag_args[[2]]), ncol = ncol(imag_args[[2]]))
  else if(is.na(imag_args[[2]]))
    imag_args[[2]] <- matrix(0, nrow = nrow(real_args[[2]]), ncol = ncol(real_args[[2]]))
  mat <- bmat(list(list(real_args[[2]], -imag_args[[2]]),
                   list(imag_args[[2]], real_args[[2]])
              ))
  return(list(copy(expr, list(vec, mat)), NA))
}

Complex2Real.param_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is_real(expr))
    return(list(expr, NA))
  else if(is_imag(expr)) {
    imag <- CallbackParam(function() { list(Im(value(expr)), dim(expr)) })
    return(list(NA, imag))
  } else {
    real <- CallbackParam(function() { list(Re(value(expr)), dim(expr)) })
    imag <- CallbackParam(function() { list(Im(value(expr)), dim(expr)) })
    return(list(real, imag))
  }
}

Complex2Real.pnorm_canon <- function(expr, real_args, imag_args, real2imag) {
  abs_args <- Complex2Real.abs_canon(expr, real_args, imag_args, real2imag)
  abs_real_args <- abs_args[[1]]
  return(list(copy(expr, list(abs_real_args)), NA))
}

Complex2Real.variable_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is_real(expr))
    return(list(expr, NA))

  imag <- Variable(dim(expr), var_id = real2imag[[as.character(id(expr))]])
  if(is_imag(expr))
    return(list(NA, imag))
  else if(is_complex(expr) && is_hermitian(expr))
    return(list(Variable(dim(expr), var_id = id(expr), symmetric = TRUE), (imag - t(imag))/2))
  else   # Complex.
    return(list(Variable(dim(expr), var_id = id(expr)), imag))
}

Complex2Real.zero_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is.na(imag_args[[1]]))
    return(list(copy(expr, real_args), NA))

  imag_cons <- ZeroConstraint(imag_args[[1]], constr_id = real2imag[[as.character(id(expr))]])
  if(is.na(real_args[[1]]))
    return(list(NA, imag_cons))
  else
    return(list(copy(expr, real_args), imag_cons))
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
                                   Zero = Complex2Real.zero_canon,
                                   PSD = Complex2Real.hermitian_canon,
                                   
                                   Abs = Complex2Real.abs_canon,
                                   Norm1 = Complex2Real.pnorm_canon,
                                   NormInf = Complex2Real.pnorm_canon,
                                   Pnorm = Complex2Real.pnorm_canon,
                                   
                                   LambdaMax = Complex2Real.hermitian_canon,
                                   LogDet = Complex2Real.norm_nuc_canon,
                                   NormNuc = Complex2Real.norm_nuc_canon,
                                   SigmaMax = Complex2Real.hermitian_canon,
                                   QuadForm = Complex2Real.quad_canon,
                                   MatrixFrac = Complex2Real.matrix_frac_canon,
                                   LambdaSumLargest = Complex2Real.lambda_sum_largest_canon)
