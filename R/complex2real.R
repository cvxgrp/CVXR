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

Complex2Real.canonicalize_expr(expr, real_args, imag_args, real2imag, leaf_map) {
  if(is(expr, "Expression") && length(variables(expr)) == 0) {
    # Parameterized expressions are evaluated in a subsequent reduction.
    if(length(parameters(expr)) > 0)
      stop("Unimplemented")
    else   # Non-parameterized expressions are evaluated immediately.
      return(elim_cplx_methods$Constant(Constant(value(expr)), real_args, imag_args, real2imag)
  } else if(type(expr) %in% names(elim_cplx_methods)) {
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
