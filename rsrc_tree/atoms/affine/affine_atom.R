## CVXPY SOURCE: cvxpy/atoms/affine/affine_atom.py
#'
#' The AffAtom class.
#'
#' This virtual class represents an affine atomic expression.
#'
#' @name AffAtom-class
#' @aliases AffAtom
#' @rdname AffAtom-class
AffAtom <- setClass("AffAtom", contains = c("VIRTUAL", "Atom"))

#' @param object An \linkS4class{AffAtom} object.
#' @describeIn AffAtom Does the atom handle complex numbers?
setMethod("allow_complex", "AffAtom", function(object) { TRUE })

#' @describeIn AffAtom The sign of the atom.
setMethod("sign_from_args", "AffAtom", function(object) { sum_signs(object@args) })

#' @describeIn AffAtom Is the atom imaginary?
setMethod("is_imag", "AffAtom", function(object) { all(sapply(object@args, is_imag)) })

#' @describeIn AffAtom Is the atom complex valued?
setMethod("is_complex", "AffAtom", function(object) { any(sapply(object@args, is_complex)) })

#' @describeIn AffAtom The atom is convex.
setMethod("is_atom_convex", "AffAtom", function(object) { TRUE })

#' @describeIn AffAtom The atom is concave.
setMethod("is_atom_concave", "AffAtom", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn AffAtom The atom is weakly increasing in every argument.
setMethod("is_incr", "AffAtom", function(object, idx) { TRUE })

#' @describeIn AffAtom The atom is not weakly decreasing in any argument.
setMethod("is_decr", "AffAtom", function(object, idx) { FALSE })

#' @describeIn AffAtom Is every argument quadratic?
setMethod("is_quadratic", "AffAtom", function(object) { all(sapply(object@args, is_quadratic)) })

#' @describeIn AffAtom Does the affine head of the expression contain a quadratic term? The affine head is all nodes with a path to the root node that does not pass through any non-affine atom. If the root node is non-affine, then the affine head is the root alone.
setMethod("has_quadratic_term", "AffAtom", function(object) { any(sapply(object@args, has_quadratic_term)) }) 

#' @describeIn AffAtom Is every argument quadratic of piecewise affine?
setMethod("is_qpwa", "AffAtom", function(object) { all(sapply(object@args, is_qpwa)) })

#' @describeIn AffAtom Is every argument piecewise linear?
setMethod("is_pwl", "AffAtom", function(object) { all(sapply(object@args, is_pwl)) })

#' @describeIn AffAtom Is the atom a positive semidefinite matrix?
setMethod("is_psd", "AffAtom", function(object) {
  for(idx in seq_len(length(object@args))) {
    arg <- object@args[[idx]]
    if(!((is_incr(object, idx) && is_psd(arg)) || (is_decr(object, idx) && is_nsd(arg))))
      return(FALSE)
  }
  return(TRUE)
})

#' @describeIn AffAtom Is the atom a negative semidefinite matrix?
setMethod("is_nsd", "AffAtom", function(object) {
  for(idx in seq_len(length(object@args))) {
    arg <- object@args[[1]]
    if(!((is_decr(object, idx) && is_psd(arg)) || (is_incr(object, idx) && is_nsd(arg))))
      return(FALSE)
  }
  return(TRUE)
})

#' @param values A list of numeric values for the arguments
#' @describeIn AffAtom Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "AffAtom", function(object, values) {
  # TODO: Should be a simple function in CVXcore for this.
  # Make a fake LinOp tree for the function
  fake_args <- list()
  var_offsets <- c()
  var_names <- c()
  offset <- 0
  for(idx in seq_len(length(object@args))) {
    arg <- object@args[[idx]]
    if(is_constant(arg))
      fake_args <- c(fake_args, list(canonical_form(Constant(value(arg)))[[1]]))
    else {
      fake_args <- c(fake_args, list(lu.create_var(dim(arg), idx)))
      var_offsets <- c(var_offsets, offset)
      var_names <- c(var_names, idx)
      offset <- offset + size(arg)
    }
  }
  var_length <- offset
  names(var_offsets) <- var_names
  graph <- graph_implementation(object, fake_args, dim(object), get_data(object))
  fake_expr <- graph[[1]]

  # Get the matrix representation of the function.
  # prob_mat <- get_problem_matrix(list(fake_expr), var_offsets)
  # V <- prob_mat[[1]]
  # I <- prob_mat[[2]] + 1   # TODO: R uses 1-indexing, but get_problem_matrix returns with 0-indexing
  # J <- prob_mat[[3]] + 1
  # dims <- c(offset, size(object))
  # stacked_grad <- sparseMatrix(i = J, j = I, x = V, dims = dims)
  
  param_to_size <- list()
  param_to_col <- list()
  param_to_size[[as.character(CONSTANT_ID)]] <- 1
  param_to_col[[as.character(CONSTANT_ID)]] <- 0
  canon_mat <- get_problem_matrix(list(fake_expr), var_length, var_offsets, param_to_size, param_to_col, size(object))
  
  # HACK TODO Convert tensors back to vectors.
  dims <- c(var_length + 1, size(object))
  stacked_grad <- matrix(t(canon_mat), nrow = dims[1], ncol = dims[2], byrow = TRUE)
  stacked_grad <- Matrix(stacked_grad, sparse = TRUE)   # TODO: How to reshape sparse matrix with byrow = TRUE?
  stacked_grad <- stacked_grad[-nrow(stacked_grad),]   # Remove last row.

  # Break up into per argument matrices.
  grad_list <- list()
  start <- 1
  for(arg in object@args) {
    if(is_constant(arg)) {
      grad_dim <- c(size(arg), dims[2])
      if(all(grad_dim == c(1,1)))
        grad_list <- c(grad_list, list(0))
      else
        grad_list <- c(grad_list, list(sparseMatrix(i = c(), j = c(), dims = grad_dim)))
    } else {
      stop <- start + size(arg)
      if(stop == start)
        grad_list <- c(grad_list, list(sparseMatrix(i = c(), j = c(), dims = c(0, dim[2]))))
      else
        grad_list <- c(grad_list, list(stacked_grad[start:(stop-1),]))
      start <- stop
    }
  }
  return(grad_list)
})
