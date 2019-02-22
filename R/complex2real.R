setClass("Complex2Real", contains = "Reduction")

setMethod("accepts", signature(object = "Complex2Real", problem = "Problem"), function(object, problem) {
  leaves <- c(variables(problem), parameters(problem), constants(problem))
  any(sapply(leaves, function(l) { is_complex(l) }))
})

setMethod("apply", signature(object = "Complex2Real", problem = "Problem"), function(object, problem) {
  inverse_data <- InverseData(problem)

  leaf_map <- list()
  obj <- canonicalize_tree(object, problem@objective, inverse_data@real2imag, leaf_map)
  real_obj <- obj[[1]]
  imag_obj <- obj[[2]]

  if(length(imag_obj) > 0)
    stop("Cannot have imaginary component in canonicalized objective")

  constrs <- list()
  for(constraint in problem@constraints) {
    constr <- canonicalize_tree(object, constraint, inverse_data@real2imag, leaf_map)
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
      } else if((is(cons, "Zero") || is(cons, "NonPos")) && is_complex(cons)) {
        imag_id <- inverse_data@real2imag[[cid]]
        dvars[[cid]] <- solution@dual_vars[[cid]] + 1i*solution@dual_vars[[imag_id]]
      # For PSD constraints.
      } else if(is(cons, "PSD") && is_complex(cons)) {
        n <- cons@args[[1]]@shape[1]
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
    outs <- canonicalize_expr(expr, real_args, imag_args, real2imag, leaf_map)
    return(list(real_out = outs[[1]], imag_out = outs[[2]]))
  }
}

Complex2Real.canonicalize_expr <- function(expr, real_args, imag_args, real2imag, leaf_map) {
  if(type(expr) %in% names(elim_cplx_methods)) {
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

abs_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is.na(real_args[[1]]))   # Imaginary
    output <- abs(imag_args[[1]])
  else if(is.na(imag_args[[1]]))   # Real
    output <- abs(real_args[[1]])
  else {   # Complex
    real <- real_args[[1]]
    imag <- imag_args[[1]]
    norms <- pnorm(vstack(c(real, imag)), p = 2, axis = 0)
    output <- reshape(norms, real_args[[1]]@shape)
  }
  return(list(output, NA))
}

# Affine canonicalization.
separable_canon <- function(expr, real_args, imag_args, real2imag) {
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

real_canon <- function(expr, real_args, imag_args, real2imag) {
  list(real_args[[1]], NA)
}

imag_canon <- function(expr, real_args, imag_args, real2imag) {
  list(imag_args[[1]], NA)
}

conj_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is.na(imag_args[[1]]))
    imag_arg <- NA
  else
    imag_arg <- -imag_args[[1]]
  return(list(real_args[[1]], imag_arg))
}

join <- function(expr, lh_arg, rh_arg) {
  # Helper function to combine arguments.
  if(is.na(lh_arg) || is.na(rh_arg))
    return(NA)
  else
    return(copy(expr, list(lh_arg, rh_arg)))
}

add <- function(lh_arg, rh_arg, neg = FALSE) {
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

binary_canon <- function(expr, real_args, imag_args, real2imag) {
  # Canonicalize functions like multiplication.
  real_by_real <- join(expr, real_args[[1]], real_args[[2]])
  imag_by_imag <- join(expr, imag_args[[1]], imag_args[[2]])
  real_by_imag <- join(expr, real_args[[1]], imag_args[[2]])
  imag_by_real <- join(expr, imag_args[[1]], real_args[[2]])
  real_output <- add(real_by_real, imag_by_imag, neg = TRUE)
  imag_output <- add(real_by_imag, imag_by_real, neg = TRUE)
  return(list(real_output, imag_output))
}

constant_canon <- function(expr, real_args, imag_args, real2imag) {
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
hermitian_canon <- function(expr, real_args, imag_args, real2imag) {
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

norm_nuc_canon <- function(expr, real_args, imag_args, real2imag) {
  # Canonicalize nuclear norm with Hermitian matrix input.
  # Divide by two because each eigenvalue is repeated twice.
  canon <- hermitian_canon(expr, real_args, imag_args, real2imag)
  real <- canon[[1]]
  imag <- canon[[2]]
  if(!is.na(imag_args[[1]]))
    real <- real/2
  return(list(real, imag))
}

lambda_sum_largest_canon <- function(expr, real_args, imag_args, real2imag) {
  # Canonicalize nuclear norm with Hermitian matrix input.
  # Divide by two because each eigenvalue is repeated twice.
  canon <- hermitian_canon(expr, real_args, imag_args, real2imag)
  real <- canon[[1]]
  imag <- canon[[2]]
  real@k <- 2*real@k
  if(!is.na(imag_args[[1]]))
    real <- real/2
  return(list(real, imag))
}

at_least_2D <- function(expr) {
  # Upcast 0D and 1D to 2D.
  if(expr@ndim < 2)
    return(reshape(expr, c(size(expr), 1)))
  else
    return(expr)
}

quad_canon <- function(expr, real_args, imag_args, real2imag) {
  # Convert quad_form to real.
  if(is.na(imag_args[[1]])) {
    vec <- real_args[[1]]
    mat <- real_args[[2]]
  } else if(is.na(real_args[[1]])) {
    vec <- imag_args[[1]]
    mat <- real_args[[2]]
  } else {
    vec <- vstack(list(at_least_2D(real_args[[1]]),
                       at_least_2D(imag_args[[1]])))
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

matrix_frac_canon <- function(expr, real_args, imag_args, real2imag) {
  # Convert matrix_frac to real.
  if(is.na(real_args[[1]]))
    real_args[[1]] <- matrix(0, nrow = nrow(imag_args[[1]]), ncol = ncol(imag_args[[1]]))
  if(is.na(imag_args[[1]]))
    imag_args[[1]] <- matrix(0, nrow = nrow(real_args[[1]]), ncol = ncol(real_args[[1]]))
  vec <- vstack(list(at_least_2D(real_args[[1]]),
                     at_least_2D(imag_args[[1]])))
  if(is.na(real_args[[2]]))
    real_args[[2]] <- matrix(0, nrow = nrow(imag_args[[2]]), ncol = ncol(imag_args[[2]]))
  else if(is.na(imag_args[[2]]))
    imag_args[[2]] <- matrix(0, nrow = nrow(real_args[[2]]), ncol = ncol(real_args[[2]]))
  mat <- bmat(list(list(real_args[[2]], -imag_args[[2]]),
                   list(imag_args[[2]], real_args[[2]])
              ))
  return(list(copy(expr, list(vec, mat)), NA))
}

param_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is_real(expr))
    return(list(expr, NA))
  else if(is_imag(expr)) {
    imag <- CallbackParam(function() { list(Im(value(expr)), shape(expr)) })
    return(list(NA, imag))
  } else {
    real <- CallbackParam(function() { list(Re(value(expr)), shape(expr)) })
    imag <- CallbackParam(function() { list(Im(value(expr)), shape(expr)) })
    return(list(real, imag))
  }
}

pnorm_canon <- function(expr, real_args, imag_args, real2imag) {
  abs_args <- abs_canon(expr, real_args, imag_args, real2imag)
  abs_real_args <- abs_args[[1]]
  return(list(copy(expr, list(abs_real_args)), NA))
}

variable_canon <- function(expr, real_args, imag_args, real2imag) {
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

zero_canon <- function(expr, real_args, imag_args, real2imag) {
  if(is.na(imag_args[[1]]))
    return(list(copy(expr, real_args), NA))

  imag_cons <- Zero(imag_args[[1]], constr_id = real2imag[[as.character(id(expr))]])
  if(is.na(real_args[[1]]))
    return(list(NA, imag_cons))
  else
    return(list(copy(expr, real_args), imag_cons))
}
