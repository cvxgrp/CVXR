## CVXPY SOURCE: cvxpy/reductions/complex2real/complex2real.py

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

Complex2Real.UNIMPLEMENTED_COMPLEX_DUALS <- c("SOC", "OpRelEntrConeQuad")

#' @param object A \linkS4class{Complex2Real} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn Complex2Real Checks whether or not the problem involves any complex numbers.
setMethod("accepts", signature(object = "Complex2Real", problem = "Problem"), function(object, problem) {
  Complex2Real.accepts(problem)
})

#' @describeIn Complex2Real Converts a Complex problem into a Real one.
setMethod("perform", signature(object = "Complex2Real", problem = "Problem"), function(object, problem) {
  inverse_data <- InverseData(problem)
  real2imag <- list()
  for(var in variables(problem)) {
    if(is_complex(var))
      real2imag[[as.character(id(var))]] <- lu.get_id()
  }
  for(cons in problem@constraints) {
    if(is_complex(cons)) {
      real2imag[[as.character(id(cons))]] <- lu.get_id()
    }
  }
  inverse_data@real2imag <- real2imag

  leaf_map <- list()
  obj <- Complex2Real.canonicalize_tree(problem@objective, inverse_data@real2imag, leaf_map)
  real_obj <- obj[[1]]
  imag_obj <- obj[[2]]

  if(length(imag_obj) > 0)
    stop("Cannot have imaginary component in canonicalized objective")

  constrs <- list()
  for(constraint in problem@constraints) {
    # real2imag maps variable id to a potential new variable created for the imaginary part.
    canon <- Complex2Real.canonicalize_tree(constraint, inverse_data@real2imag, leaf_map)
    real_constrs <- canon[[1]]
    imag_constrs <- canon[[2]]
    if(is.list(real_constrs) || is(real_constrs, "Constraint"))
      constrs <- c(constrs, real_constrs)
    if(is.list(imag_constrs) || is(imag_constrs, "Constraint"))
      constrs <- c(constrs, imag_constrs)
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
    #
    # Primal variables
    #
    for(vid in names(inverse_data@id2var)) {
      var <- inverse_data@id2var[[vid]]
      if(is_real(var))   # Purely real variables
        pvars[[vid]] <- solution@primal_vars[[vid]]
      else if(is_imag(var)) {   # Purely imaginary variables
        imag_id <- inverse_data@real2imag[[vid]]
        pvars[[vid]] <- 1i*solution@primal_vars[[as.character(imag_id)]]
      } else if(is_complex(var) && is_hermitian(var)) {   # Hermitian variables
        pvars[[vid]] <- solution@primal_vars[[vid]]
        imag_id <- inverse_data@real2imag[[vid]]

        if(as.character(imag_id) %in% names(solution@primal_vars)) {
          imag_val <- solution@primal_vars[[as.character(imag_id)]]
          imag_val <- value(UpperTri.vec_to_upper_tri(imag_val, TRUE))
          imag_val <- imag_val - t(imag_val)
          pvars[[vid]] <- pvars[[vid]] + 1i*imag_val
        }
      } else if(is_complex(var)) {   # General complex variables
        pvars[[vid]] <- solution@primal_vars[[vid]]
        imag_id <- inverse_data@real2imag[[vid]]
        if(as.character(imag_id) %in% names(solution@primal_vars)) {
          imag_val <- solution@primal_vars[[as.character(imag_id)]]
          pvars[[vid]] <- pvars[[vid]] + 1i*imag_val
        }
      }
    }

    if(!is.null(solution@dual_vars)) {
      #
      # Dual variables
      #
      for(cid in names(inverse_data@id2cons)) {
        cons <- inverse_data@id2cons[[cid]]
        if(is_real(cons))
          dvars[[cid]] <- solution@dual_vars[[cid]]
        else if(is_imag(cons)) {
          imag_id <- inverse_data@real2imag[[cid]]
          dvars[[cid]] <- 1i*solution@dual_vars[[as.character(imag_id)]]
        # All cases that follow are for complex-valued constraints:
        #     1. Check inequality / equality constraints.
        #     2. Check PSD constraints.
        #     3. Check if a constraint is known to lack a complex dual implementation.
        #     $. Raise an error.
        } else if(is(cons, "ZeroConstraint") || is(cons, "EqConstraint") || is(cons, "IneqConstraint") || is(cons, "NonNegConstraint") || is(cons, "NonPosConstraint")) {
          imag_id <- inverse_data@real2imag[[cid]]
          if(imag_id %in% names(solution@dual_vars))
            dvars[[cid]] <- solution@dual_vars[[cid]] + 1i*solution@dual_vars[[imag_id]]
          else
            dvars[[cid]] <- solution@dual_vars[[cid]]
        } else if(is(cons, "PSDConstraint")) {
          # Suppose we have a constraint con_x = X >> 0 where X is Hermitian.
          #
          # Define the matrix
          #     Y := [[ re(X), im(X)],
          #           [-im(X), re(X)]]
          # and the constraint con_y = Y >> 0.
          #
          # The real part of the dual variable for con_x is the upper-left block of the dual variable for con_y.
          # The imaginary part of the dual variable for con_x is the upper-right block of the dual variable for con_y.
          n <- nrow(cons@args[[1]])
          dual <- solution@dual_vars[[cid]]
          dvars[[cid]] <- dual[1:n,1:n] + 1i*dual[(n+1):nrow(dual),1:n]
        } else if(any(sapply(Complex2Real.UNIMPLEMENTED_COMPLEX_DUALS, function(c) { isinstance(cons, c) } ))) {
          # TODO: Implement dual variable recovery
        } else
          stop("Unknown constraint type")
      }
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
  if(inherits(expr, names(Complex2Real.CANON_METHODS))) {
    expr_id <- as.character(expr@id)
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
