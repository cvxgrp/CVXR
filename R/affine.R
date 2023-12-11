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
      fake_args <- c(fake_args, list(create_var(dim(arg), idx)))
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

#'
#' The AddExpression class.
#'
#' This class represents the sum of any number of expressions.
#'
#' @slot arg_groups A \code{list} of \linkS4class{Expression}s and numeric data.frame, matrix, or vector objects.
#' @name AddExpression-class
#' @aliases AddExpression
#' @rdname AddExpression-class
.AddExpression <- setClass("AddExpression", representation(arg_groups = "list"), prototype(arg_groups = list()), contains = "AffAtom")
AddExpression <- function(arg_groups = list()) { .AddExpression(arg_groups = arg_groups) }

setMethod("initialize", "AddExpression", function(.Object, ..., arg_groups = list()) {
  .Object@arg_groups <- arg_groups
  .Object <- callNextMethod(.Object, ..., atom_args = arg_groups)   # Casts R values to Constant objects
  .Object@args <- lapply(.Object@args, function(group) { if(is(group,"AddExpression")) group@args else group })
  .Object@args <- flatten_list(.Object@args)   # Need to flatten list of expressions
  .Object
})

#' @param x,object An \linkS4class{AddExpression} object.
#' @describeIn AddExpression The dimensions of the expression.
setMethod("dim_from_args", "AddExpression", function(object) { sum_dims(lapply(object@args, dim)) })

setMethod("expand_args", "AddExpression", function(expr) { 
  if(is(expr, "AddExpression"))
    return(expr@args)
  else
    return(list(expr))
})

#' @describeIn AddExpression The string form of the expression.
setMethod("name", "AddExpression", function(x) {
  paste(sapply(x@args, name), collapse = " + ")
})

#' @param values A list of arguments to the atom.
#' @describeIn AddExpression Sum all the values.
setMethod("to_numeric", "AddExpression", function(object, values) {
  values <- lapply(values, intf_convert_if_scalar)
  Reduce("+", values)
})

#' @describeIn AddExpression Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "AddExpression", function(object) { TRUE })

#' @describeIn AddExpression Is the atom log-log convex?
setMethod("is_atom_log_log_concave", "AddExpression", function(object) { FALSE })

#' @describeIn AddExpression Is the atom symmetric?
setMethod("is_symmetric", "AddExpression", function(object) {
  symm_args <- all(sapply(object@args, is_symmetric))
  return(dim(object)[1] == dim(object)[2] && symm_args)
})

#' @describeIn AddExpression Is the atom hermitian?
setMethod("is_hermitian", "AddExpression", function(object) {
  herm_args <- all(sapply(object@args, is_hermitian))
  return(dim(object)[1] == dim(object)[2] && herm_args)
})

# As initialize takes in the arg_groups instead of args, we need a special copy function.
#' @param args An optional list of arguments to reconstruct the atom. Default is to use current args of the atom.
#' @param id_objects Currently unused.
#' @describeIn AddExpression Returns a shallow copy of the AddExpression atom
setMethod("copy", "AddExpression", function(object, args = NULL, id_objects = list()) {
  if(is.null(args))
    args <- object@arg_groups
  do.call(class(object), list(arg_groups = args))
})

AddExpression.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  arg_objs <- lapply(arg_objs, function(arg) {
    if(!all(arg$dim == dim) && lo.is_scalar(arg)) 
      lo.promote(arg, dim) 
    else 
      arg })
  list(lo.sum_expr(arg_objs), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn AddExpression The graph implementation of the expression.
setMethod("graph_implementation", "AddExpression", function(object, arg_objs, dim, data = NA_real_) {
  AddExpression.graph_implementation(arg_objs, dim, data)
})

#'
#' The UnaryOperator class.
#'
#' This base class represents expressions involving unary operators.
#'
#' @slot expr The \linkS4class{Expression} that is being operated upon.
#' @slot op_name A \code{character} string indicating the unary operation.
#' @name UnaryOperator-class
#' @aliases UnaryOperator
#' @rdname UnaryOperator-class
UnaryOperator <- setClass("UnaryOperator", representation(expr = "Expression"), contains = "AffAtom")

setMethod("initialize", "UnaryOperator", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

setMethod("op_name", "UnaryOperator", function(object) { stop("Unimplemented") })
setMethod("op_func", "UnaryOperator", function(object) { stop("Unimplemented") })

#' @param x,object A \linkS4class{UnaryOperator} object.
#' @describeIn UnaryOperator Returns the expression in string form.
setMethod("name", "UnaryOperator", function(x) {
  paste(op_name(x), name(x@args[[1]]), sep = "")
})

#' @param values A list of arguments to the atom.
#' @describeIn UnaryOperator Applies the unary operator to the value.
setMethod("to_numeric", "UnaryOperator", function(object, values) {
  op_func(object)(values[[1]])
})

#'
#' The NegExpression class.
#'
#' This class represents the negation of an affine expression.
#'
#' @name NegExpression-class
#' @aliases NegExpression
#' @rdname NegExpression-class
.NegExpression <- setClass("NegExpression", contains = "UnaryOperator")
NegExpression <- function(expr) { .NegExpression(expr = expr) }

setMethod("op_name", "NegExpression", function(object) { "-" })
setMethod("op_func", "NegExpression", function(object) { function(x) { -x } })

#' @param object A \linkS4class{NegExpression} object.
#' @describeIn NegExpression The (row, col) dimensions of the expression.
setMethod("dim_from_args", "NegExpression", function(object) { dim(object@args[[1]]) })

#' @describeIn NegExpression The (is positive, is negative) sign of the expression.
setMethod("sign_from_args", "NegExpression", function(object) { c(is_nonpos(object@args[[1]]), is_nonneg(object@args[[1]])) })

#' @param idx An index into the atom.
#' @describeIn NegExpression The expression is not weakly increasing in any argument.
setMethod("is_incr", "NegExpression", function(object, idx) { FALSE })

#' @describeIn NegExpression The expression is weakly decreasing in every argument.
setMethod("is_decr", "NegExpression", function(object, idx) { TRUE })

#' @describeIn NegExpression Is the expression symmetric?
setMethod("is_symmetric", "NegExpression", function(object) { is_symmetric(object@args[[1]]) })

#' @describeIn NegExpression Is the expression Hermitian?
setMethod("is_hermitian", "NegExpression", function(object) { is_hermitian(object@args[[1]]) })

NegExpression.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.neg_expr(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn NegExpression The graph implementation of the expression.
setMethod("graph_implementation", "NegExpression", function(object, arg_objs, dim, data = NA_real_) {
  NegExpression.graph_implementation(arg_objs, dim, data)
})

#' The BinaryOperator class.
#'
#' This base class represents expressions involving binary operators.
#'
#' @slot lh_exp The \linkS4class{Expression} on the left-hand side of the operator.
#' @slot rh_exp The \linkS4class{Expression} on the right-hand side of the operator.
#' @slot op_name A \code{character} string indicating the binary operation.
#' @name BinaryOperator-class
#' @aliases BinaryOperator
#' @rdname BinaryOperator-class
BinaryOperator <- setClass("BinaryOperator", representation(lh_exp = "ConstValORExpr", rh_exp = "ConstValORExpr"), contains = "AffAtom")

setMethod("initialize", "BinaryOperator", function(.Object, ..., lh_exp, rh_exp) {
  .Object@lh_exp = lh_exp
  .Object@rh_exp = rh_exp
  callNextMethod(.Object, ..., atom_args = list(.Object@lh_exp, .Object@rh_exp))
})

setMethod("op_name", "BinaryOperator", function(object) { "BINARY_OP" })

#' @param x,object A \linkS4class{BinaryOperator} object.
#' @describeIn BinaryOperator Returns the name of the BinaryOperator object.
setMethod("name", "BinaryOperator", function(x) {
  pretty_args <- list()
  for(a in x@args) {
    if(is(a, "AddExpression") || is(a, "DivExpression"))
      pretty_args <- list(pretty_args, paste("(", name(a), ")", sep = ""))
    else
      pretty_args <- list(pretty_args, name(a))
  }
  paste(pretty_args[[1]], op_name(x), pretty_args[[2]])
})

#' @param values A list of arguments to the atom.
#' @describeIn BinaryOperator Apply the binary operator to the values.
setMethod("to_numeric", "BinaryOperator", function(object, values) {
  values <- lapply(values, intf_convert_if_scalar)
  Reduce(op_name(object), values)
})

#' @describeIn BinaryOperator Default to rule for multiplication.
setMethod("sign_from_args", "BinaryOperator", function(object) { mul_sign(object@args[[1]], object@args[[2]]) })

#' @describeIn BinaryOperator Is the expression imaginary?
setMethod("is_imag", "BinaryOperator", function(object) {
  (is_imag(object@args[[1]]) && is_real(object@args[[2]])) || (is_real(object@args[[1]]) && is_imag(object@args[[2]]))
})

#' @describeIn BinaryOperator Is the expression complex valued?
setMethod("is_complex", "BinaryOperator", function(object) {
  (is_complex(object@args[[1]]) || is_complex(object@args[[2]])) && !(is_imag(object@args[[1]]) && is_imag(object@args[[2]]))
})

#'
#' The MulExpression class.
#'
#' This class represents the matrix product of two linear expressions.
#' See \linkS4class{Multiply} for the elementwise product.
#'
#' @seealso \linkS4class{Multiply}
#' @name MulExpression-class
#' @aliases MulExpression
#' @rdname MulExpression-class
.MulExpression <- setClass("MulExpression", contains = "BinaryOperator")
MulExpression <- function(lh_exp, rh_exp) { .MulExpression(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("op_name", "MulExpression", function(object) { "*" })

#' @param object A \linkS4class{MulExpression} object.
#' @param values A list of arguments to the atom.
#' @describeIn MulExpression Matrix multiplication.
setMethod("to_numeric", "MulExpression", function(object, values) {
  if(is.null(dim(object@args[[1]])) || is.null(dim(object@args[[2]])))
    return(values[[1]] * values[[2]])
  else
    return(values[[1]] %*% values[[2]])
})

#' @describeIn MulExpression The (row, col) dimensions of the expression.
setMethod("dim_from_args", "MulExpression", function(object) { mul_dims(dim(object@args[[1]]), dim(object@args[[2]])) })

#' @describeIn MulExpression Multiplication is convex (affine) in its arguments only if one of the arguments is constant.
setMethod("is_atom_convex", "MulExpression", function(object) {
  if(dpp_scope_active()) {
    # This branch applies curvature rules for DPP.
    #
    # Because a DPP scope is active, parameters will be
    # treated as affine (like variables, not constants) by curvature
    # analysis methods.
    #
    # Like under DCP, a product x * y is convex if x or y is constant.
    # If neither x nor y is constant, then the product is DPP
    # if one of the expressions is affine in its parameters and the
    # other is parameter-free.
    x <- object@args[[1]]
    y <- object@args[[2]]
    (is_constant(x) || is_constant(y)) || (is_param_affine(x) && is_param_free(y)) || (is_param_affine(y) && is_param_free(x))
  } else
    is_constant(object@args[[1]]) || is_constant(object@args[[2]])
})

#' @describeIn MulExpression If the multiplication atom is convex, then it is affine.
setMethod("is_atom_concave", "MulExpression", function(object) { is_atom_convex(object) })

#' @describeIn MulExpression Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "MulExpression", function(object) { TRUE })

#' @describeIn MulExpression Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "MulExpression", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn MulExpression Is the left-hand expression positive?
setMethod("is_incr", "MulExpression", function(object, idx) { is_nonneg(object@args[[3-idx]]) })

#' @describeIn MulExpression Is the left-hand expression negative?
setMethod("is_decr", "MulExpression", function(object, idx) { is_nonpos(object@args[[3-idx]]) })

#' @param values A list of numeric values for the arguments
#' @describeIn MulExpression Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "MulExpression", function(object, values) {
  if(is_constant(object@args[[1]]) || is_constant(object@args[[2]]))
    return(callNextMethod(object, values))
  
  # TODO: Verify that the following code is correct for non-affine arguments.
  X <- values[[1]]
  Y <- values[[2]]
  
  DX_rows <- size(object@args[[1]])
  block_rows <- DX_rows/nrow(Y)
  
  # DX = [diag(Y11), diag(Y12), ...]
  #      [diag(Y21), diag(Y22), ...]
  #      [   ...       ...      ...]
  DX <- kronecker(Y, sparseMatrix(i = 1:block_rows, j = 1:block_rows, x = 1))
  cols <- ifelse(length(dim(object@args[[2]])) == 1, 1, ncol(object@args[[2]]))
  DY <- Matrix(bdiag(lapply(1:cols, function(k) { t(X) })), sparse = TRUE)
  return(list(DX, DY))
})

MulExpression.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  # Promote the right-hand side to a diagonal matrix if necessary
  lhs <- arg_objs[[1]]
  rhs <- arg_objs[[2]]
  if(lo.is_const(lhs))
    return(list(lo.mul_expr(lhs, rhs, dim), list()))
  else if(lo.is_const(rhs))
    return(list(lo.rmul_expr(lhs, rhs, dim), list()))
  else
    stop("Product of two non-constant expressions is not DCP.")
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn MulExpression The graph implementation of the expression.
setMethod("graph_implementation", "MulExpression", function(object, arg_objs, dim, data = NA_real_) {
  MulExpression.graph_implementation(arg_objs, dim, data)
})

#'
#' The Multiply class.
#'
#' This class represents the elementwise product of two expressions.
#'
#' @name Multiply-class
#' @aliases Multiply
#' @rdname Multiply-class
.Multiply <- setClass("Multiply", contains = "MulExpression")

#' @param lh_exp An \linkS4class{Expression} or R numeric data.
#' @param rh_exp An \linkS4class{Expression} or R numeric data.
#' @rdname Multiply-class
Multiply <- function(lh_exp, rh_exp) { .Multiply(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("initialize", "Multiply", function(.Object, ..., lh_exp, rh_exp) {
  lh_exp <- as.Constant(lh_exp)
  rh_exp <- as.Constant(rh_exp)
  if(is_scalar(lh_exp) && !is_scalar(rh_exp))
    lh_exp <- promote(lh_exp, dim(rh_exp))
  else if(is_scalar(rh_exp) && !is_scalar(lh_exp))
    rh_exp <- promote(rh_exp, dim(lh_exp))
  callNextMethod(.Object, ..., lh_exp = lh_exp, rh_exp = rh_exp)
})

#' @param object A \linkS4class{Multiply} object.
#' @param values A list of arguments to the atom.
#' @describeIn Multiply Multiplies the values elementwise.
setMethod("to_numeric", "Multiply", function(object, values) { values[[1]] * values[[2]] })

#' @describeIn Multiply The sum of the argument dimensions - 1.
setMethod("dim_from_args", "Multiply", function(object) { sum_dims(lapply(object@args, dim)) })

#' @describeIn Multiply Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Multiply", function(object) { TRUE })

#' @describeIn Multiply Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Multiply", function(object) { TRUE })

#' @describeIn Multiply Is the atom quasiconvex?
setMethod("is_atom_quasiconvex", "Multiply", function(object) {
  (is_constant(object@args[[1]]) || is_constant(object@args[[2]])) ||
    (is_nonneg(object@args[[1]]) && is_nonpos(object@args[[2]])) ||
    (is_nonpos(object@args[[1]]) && is_nonneg(object@args[[2]]))
})

#' @describeIn Multiply Is the atom quasiconcave?
setMethod("is_atom_quasiconcave", "Multiply", function(object) {
  is_constant(object@args[[1]]) || is_constant(object@args[[2]]) ||
    all(sapply(object@args, is_nonneg)) || all(sapply(object@args, is_nonpos))
})

#' @describeIn Multiply Is the expression a positive semidefinite matrix?
setMethod("is_psd", "Multiply", function(object) {
  (is_psd(object@args[[1]]) && is_psd(object@args[[2]])) || (is_nsd(object@args[[1]]) && is_nsd(object@args[[2]]))
})

#' @describeIn Multiply Is the expression a negative semidefinite matrix?
setMethod("is_nsd", "Multiply", function(object) {
  (is_psd(object@args[[1]]) && is_nsd(object@args[[2]])) || (is_nsd(object@args[[1]]) && is_psd(object@args[[2]]))
})

Multiply.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  lhs <- arg_objs[[1]]
  rhs <- arg_objs[[2]]
  if(lo.is_const(lhs))
    return(list(lo.multiply(lhs, rhs), list()))
  else if(lo.is_const(rhs))
    return(list(lo.multiply(rhs, lhs), list()))
  else
    stop("Product of two non-constant expressions is not DCP.")
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Multiply The graph implementation of the expression.
setMethod("graph_implementation", "Multiply", function(object, arg_objs, dim, data = NA_real_) {
  Multiply.graph_implementation(arg_objs, dim, data)
})

#'
#' The DivExpression class.
#'
#' This class represents one expression divided by another expression.
#'
#' @name DivExpression-class
#' @aliases DivExpression
#' @rdname DivExpression-class
.DivExpression <- setClass("DivExpression", contains = "Multiply")
DivExpression <- function(lh_exp, rh_exp) { .DivExpression(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("op_name", "DivExpression", function(object) { "/" })

#' @param object A \linkS4class{DivExpression} object.
#' @param values A list of arguments to the atom.
#' @describeIn DivExpression Matrix division by a scalar.
setMethod("to_numeric", "DivExpression", function(object, values) {
  if(!is.null(dim(values[[2]])) && prod(dim(values[[2]])) == 1)
    values[[2]] <- as.vector(values[[2]])[1]
  return(values[[1]] / values[[2]])
})

#' @param object A \linkS4class{DivExpression} object.
#' @describeIn DivExpression Is the left-hand expression quadratic and the right-hand expression constant?
setMethod("is_quadratic", "DivExpression", function(object) {
  is_quadratic(object@args[[1]]) && is_constant(object@args[[2]])
})

#' @describeIn DivExpression Can be a quadratic term if divisor is constant.
setMethod("has_quadratic_term", "DivExpression", function(object) {
  has_quadratic_term(object@args[[1]]) && is_constant(object@args[[2]])
})

#' @describeIn DivExpression Is the expression quadratic of piecewise affine?
setMethod("is_qpwa", "DivExpression", function(object) {
  is_qpwa(object@args[[1]]) && is_constant(object@args[[2]])
})

#' @describeIn DivExpression The (row, col) dimensions of the left-hand expression.
setMethod("dim_from_args", "DivExpression", function(object) { dim(object@args[[1]]) })

#' @describeIn DivExpression Division is convex (affine) in its arguments only if the denominator is constant.
setMethod("is_atom_convex", "DivExpression", function(object) { is_constant(object@args[[2]]) })

#' @describeIn DivExpression Division is concave (affine) in its arguments only if the denominator is constant.
setMethod("is_atom_concave", "DivExpression", function(object) { is_atom_convex(object) })

#' @describeIn DivExpression Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "DivExpression", function(object) { TRUE })

#' @describeIn DivExpression Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "DivExpression", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn DivExpression Is the right-hand expression positive?
setMethod("is_incr", "DivExpression", function(object, idx) {
  if(idx == 1)
    return(is_nonneg(object@args[[2]]))
  else
    return(is_nonpos(object@args[[1]]))
})

#' @describeIn DivExpression Is the right-hand expression negative?
setMethod("is_decr", "DivExpression", function(object, idx) {
  if(idx == 1)
    return(is_nonpos(object@args[[2]]))
  else
    return(is_nonneg(object@args[[1]]))
})

DivExpression.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.div_expr(arg_objs[[1]], arg_objs[[2]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn DivExpression The graph implementation of the expression.
setMethod("graph_implementation", "DivExpression", function(object, arg_objs, dim, data = NA_real_) {
  DivExpression.graph_implementation(arg_objs, dim, data)
})

#' @param x An \linkS4class{Expression} or R numeric data.
#' @param y An \linkS4class{Expression} or R numeric data.
#' @describeIn BinaryOperator An \linkS4class{Expression} representing the standard inner product (or "scalar" product) of (x,y), conjugate-linear in x.
scalar_product <- function(x, y) {
  x <- Reshape.deep_flatten(x)
  y <- Reshape.deep_flatten(y)
  prod <- Multiply(Conjugate(x), y)
  return(SumEntries(prod))
}

#'
#' The Conjugate class.
#' 
#' This class represents the complex conjugate of an expression.
#' 
#' @slot expr An \linkS4class{Expression} or R numeric data.
#' @name Conjugate-class
#' @aliases Conjugate
#' @rdname Conjugate-class
.Conjugate <- setClass("Conjugate", representation(expr = "ConstValORExpr"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} or R numeric data.
#' @rdname Conjugate-class
Conjugate <- function(expr) { .Conjugate(expr = expr) }

setMethod("initialize", "Conjugate", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{Conjugate} object.
#' @param values A list of arguments to the atom.
#' @describeIn Conjugate Elementwise complex conjugate of the constant.
setMethod("to_numeric", "Conjugate", function(object, values) { Conj(values[[1]]) })

#' @describeIn Conjugate The (row, col) dimensions of the expression.
setMethod("dim_from_args", "Conjugate", function(object) { dim(object@args[[1]]) })

#' @param idx An index into the atom.
#' @describeIn Conjugate Is the composition weakly increasing in argument idx?
setMethod("is_incr", "Conjugate", function(object, idx) { FALSE })

#' @describeIn Conjugate Is the composition weakly decreasing in argument idx?
setMethod("is_decr", "Conjugate", function(object, idx) { FALSE })

#' @describeIn Conjugate Is the expression symmetric?
setMethod("is_symmetric", "Conjugate", function(object) { is_symmetric(object@args[[1]]) })

#' @describeIn Conjugate Is the expression hermitian?
setMethod("is_hermitian", "Conjugate", function(object) { is_hermitian(object@args[[1]]) })

Conjugate.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(arg_objs[[1]], list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Conv The graph implementation of the atom.
setMethod("graph_implementation", "Conjugate", function(object, arg_objs, dim, data = NA_real_) {
  Conjugate.graph_implementation(arg_objs, dim, data)
})

#'
#' The Conv class.
#'
#' This class represents the 1-D discrete convolution of two vectors.
#'
#' @slot lh_exp An \linkS4class{Expression} or R numeric data representing the left-hand vector.
#' @slot rh_exp An \linkS4class{Expression} or R numeric data representing the right-hand vector.
#' @name Conv-class
#' @aliases Conv
#' @rdname Conv-class
.Conv <- setClass("Conv", representation(lh_exp = "ConstValORExpr", rh_exp = "ConstValORExpr"), contains = "AffAtom")

#' @param lh_exp An \linkS4class{Expression} or R numeric data representing the left-hand vector.
#' @param rh_exp An \linkS4class{Expression} or R numeric data representing the right-hand vector.
#' @rdname Conv-class
Conv <- function(lh_exp, rh_exp) { .Conv(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("initialize", "Conv", function(.Object, ..., lh_exp, rh_exp) {
  .Object@lh_exp <- lh_exp
  .Object@rh_exp <- rh_exp
  callNextMethod(.Object, ..., atom_args = list(.Object@lh_exp, .Object@rh_exp))
})

#' @param object A \linkS4class{Conv} object.
#' @param values A list of arguments to the atom.
#' @describeIn Conv The convolution of the two values.
setMethod("to_numeric", "Conv", function(object, values) {
  .Call('_CVXR_cpp_convolve', PACKAGE = 'CVXR', as.vector(values[[1]]), as.vector(values[[2]]))
})

#' @describeIn Conv Check both arguments are vectors and the first is a constant.
setMethod("validate_args", "Conv", function(object) {
  if(!is_vector(object@args[[1]]) || !is_vector(object@args[[2]]))
    stop("The arguments to Conv must resolve to vectors.")
  if(!is_constant(object@args[[1]]))
    stop("The first argument to Conv must be constant.")
})

#' @describeIn Conv The dimensions of the atom.
setMethod("dim_from_args", "Conv", function(object) {
  lh_length <- dim(object@args[[1]])[1]
  rh_length <- dim(object@args[[2]])[1]
  c(lh_length + rh_length - 1, 1)
})

#' @describeIn Conv The sign of the atom.
setMethod("sign_from_args", "Conv", function(object) { mul_sign(object@args[[1]], object@args[[2]]) })

#' @param idx An index into the atom.
#' @describeIn Conv Is the left-hand expression positive?
setMethod("is_incr", "Conv", function(object, idx) { is_nonneg(object@args[[1]]) })

#' @param idx An index into the atom.
#' @describeIn Conv Is the left-hand expression negative?
setMethod("is_decr", "Conv", function(object, idx) { is_nonpos(object@args[[1]]) })

Conv.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.conv(arg_objs[[1]], arg_objs[[2]], dim), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Conv The graph implementation of the atom.
setMethod("graph_implementation", "Conv", function(object, arg_objs, dim, data = NA_real_) {
  Conv.graph_implementation(arg_objs, dim, data)
})

#
# Difference Matrix
#
# Returns a sparse matrix representation of the first order difference operator.
#
# @param dim The length of the matrix dimensions.
# @param axis The axis to take the difference along.
# @return A square matrix representing the first order difference.
get_diff_mat <- function(dim, axis) {
  # Construct a sparse matrix representation
  val_arr <- c()
  row_arr <- c()
  col_arr <- c()
  
  for(i in 1:dim) {
    val_arr <- c(val_arr, 1)
    row_arr <- c(row_arr, i)
    col_arr <- c(col_arr, i)
    if(i > 1) {
      val_arr <- c(val_arr, -1)
      row_arr <- c(row_arr, i)
      col_arr <- c(col_arr, i-1)
    }
  }
  
  mat <- sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = c(dim, dim))
  
  if(axis == 2)
    mat
  else
    t(mat)
}

#'
#' The CumSum class.
#'
#' This class represents the cumulative sum.
#'
#' @slot expr An \linkS4class{Expression} to be summed.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, and \code{2} indicates columns. The default is \code{2}.
#' @name CumSum-class
#' @aliases CumSum
#' @rdname CumSum-class
.CumSum <- setClass("CumSum", prototype = prototype(axis = 2), contains = c("AffAtom", "AxisAtom"))

#' @param expr An \linkS4class{Expression} to be summed.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, and \code{2} indicates columns. The default is \code{2}.
#' @rdname CumSum-class
CumSum <- function(expr, axis = 2) { .CumSum(expr = expr, axis = axis) }

#' @param object A \linkS4class{CumSum} object.
#' @param values A list of arguments to the atom.
#' @describeIn CumSum The cumulative sum of the values along the specified axis.
setMethod("to_numeric", "CumSum", function(object, values) {
  # apply(values[[1]], object@axis, base::cumsum)
  if(object@axis == 1)
    do.call(rbind, lapply(seq_len(nrow(values[[1]])), function(i) { base::cumsum(values[[1]][i,]) }))
  else if(object@axis == 2)
    do.call(cbind, lapply(seq_len(ncol(values[[1]])), function(j) { base::cumsum(values[[1]][,j]) }))
  else
    base::cumsum(values[[1]])
})

#' @describeIn CumSum The dimensions of the atom.
setMethod("dim_from_args", "CumSum", function(object) { dim(object@args[[1]]) })

#' @describeIn CumSum Returns the axis along which the cumulative sum is taken.
setMethod("get_data", "CumSum", function(object) { list(object@axis) })

#' @param values A list of numeric values for the arguments
#' @describeIn CumSum Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "CumSum", function(object, values) {
  # TODO: This is inefficient
  val_dim <- dim(values[[1]])
  collapse <- setdiff(1:length(val_dim), object@axis)
  dim <- val_dim[collapse]
  mat <- matrix(0, nrow = dim, ncol = dim)
  mat[lower.tri(mat, diag = TRUE)] <- 1

  # var <- Variable(dim(object@args[[1]]))
  var <- new("Variable", dim = dim(object@args[[1]]))
  if(object@axis == 2)
    grad <- .grad(new("MulExpression", lh_exp = mat, rh_exp = var), values)[[2]]
  else
    grad <- .grad(new("MulExpression", lh_exp = var, rh_exp = t(mat)), values)[[1]]
  list(grad)
})

CumSum.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  # Implicit O(n) definition:
  # X = Y[:1,:] - Y[1:,:]
  Y <- create_var(dim)
  axis <- data[[1]]
  collapse <- setdiff(1:length(dim), axis)
  new_dim <- dim[collapse]
  diff_mat <- get_diff_mat(new_dim, axis)
  diff_mat <- create_const(diff_mat, c(new_dim, new_dim), sparse = TRUE)

  if(axis == 2)
    diff <- lo.mul_expr(diff_mat, Y)
  else
    diff <- lo.rmul_expr(Y, diff_mat)
  list(Y, list(create_eq(arg_objs[[1]], diff)))
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn CumSum The graph implementation of the atom.
setMethod("graph_implementation", "CumSum", function(object, arg_objs, dim, data = NA_real_) {
  CumSum.graph_implementation(arg_objs, dim, data)
})

#'
#' The DiagVec class.
#'
#' This class represents the conversion of a vector into a diagonal matrix.
#'
#' @slot expr An \linkS4class{Expression} representing the vector to convert.
#' @name DiagVec-class
#' @aliases DiagVec
#' @rdname DiagVec-class
.DiagVec <- setClass("DiagVec", representation(expr = "Expression"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} representing the vector to convert.
#' @rdname DiagVec-class
DiagVec <- function(expr) { .DiagVec(expr = expr) }

setMethod("initialize", "DiagVec", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{DiagVec} object.
#' @param values A list of arguments to the atom.
#' @describeIn DiagVec Convert the vector constant into a diagonal matrix.
setMethod("to_numeric", "DiagVec", function(object, values) { diag(as.vector(values[[1]])) })

#' @describeIn DiagVec The dimensions of the atom.
setMethod("dim_from_args", "DiagVec", function(object) {
  rows <- dim(object@args[[1]])[1]
  c(rows, rows)
})

#' @describeIn DiagVec Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "DiagVec", function(object) { TRUE })

#' @describeIn DiagVec Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "DiagVec", function(object) { TRUE })

#' @describeIn DiagVec Is the expression symmetric?
setMethod("is_symmetric", "DiagVec", function(object) { TRUE })

#' @describeIn DiagVec Is the expression hermitian?
setMethod("is_hermitian", "DiagVec", function(object) { TRUE })

#' @describeIn DiagVec Is the atom positive semidefinite?
setMethod("is_psd", "DiagVec", function(object) { is_nonneg(object) })

#' @describeIn DiagVec Is the atom negative semidefinite?
setMethod("is_nsd", "DiagVec", function(object) { is_nonpos(object) })

DiagVec.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.diag_vec(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn DiagVec The graph implementation of the atom.
setMethod("graph_implementation", "DiagVec", function(object, arg_objs, dim, data = NA_real_) {
  DiagVec.graph_implementation(arg_objs, dim, data)
})

#'
#' The DiagMat class.
#'
#' This class represents the extraction of the diagonal from a square matrix.
#'
#' @slot expr An \linkS4class{Expression} representing the matrix whose diagonal we are interested in.
#' @name DiagMat-class
#' @aliases DiagMat
#' @rdname DiagMat-class
.DiagMat <- setClass("DiagMat", representation(expr = "Expression"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} representing the matrix whose diagonal we are interested in.
#' @rdname DiagMat-class
DiagMat <- function(expr) { .DiagMat(expr = expr) }

setMethod("initialize", "DiagMat", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{DiagMat} object.
#' @param values A list of arguments to the atom.
#' @describeIn DiagMat Extract the diagonal from a square matrix constant.
setMethod("to_numeric", "DiagMat", function(object, values) { diag(values[[1]]) })

#' @describeIn DiagMat The size of the atom.
setMethod("dim_from_args", "DiagMat", function(object) {
  rows <- dim(object@args[[1]])[1]
  c(rows, 1)
})

#' @describeIn DiagMat Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "DiagMat", function(object) { TRUE })

#' @describeIn DiagMat Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "DiagMat", function(object) { TRUE })

#' @describeIn DiagMat A logical value indicating whether the atom is nonnegative.
setMethod("is_nonneg", "DiagMat", function(object) {
  is_nonneg(object@args[[1]]) || is_psd(object@args[[1]])
})

DiagMat.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.diag_mat(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn DiagMat The graph implementation of the atom.
setMethod("graph_implementation", "DiagMat", function(object, arg_objs, dim, data = NA_real_) {
  DiagMat.graph_implementation(arg_objs, dim, data)
})

#' 
#' The Diag atom.
#' 
#' Turns an expression into a DiagVec object
#' 
#' @param expr An \linkS4class{Expression} that represents a vector or square matrix.
#' @return An \linkS4class{Expression} representing the diagonal vector/matrix.
#' @rdname Diag-int
Diag <- function(expr) {
  expr <- as.Constant(expr)
  if(is_vector(expr))
    return(DiagVec(Vec(expr)))
  else if(ndim(expr) == 2 && nrow(expr) == ncol(expr))
    return(DiagMat(expr = expr))
  else
    stop("Argument to Diag must be a vector or square matrix.")
}

#' 
#' The Diff atom.
#' 
#' Takes the k-th order differences
#' 
#' @param lag The degree of lag between differences
#' @param k The integer value of the order of differences
#' @param x An \linkS4class{Expression} that represents a vector
#' @param axis The axis along which to apply the function. For a 2D matrix, \code{1} indicates rows and \code{2} indicates columns.
#' @return Takes in a vector of length n and returns a vector of length n-k of the kth order differences
#' @rdname Diff-int
Diff <- function(x, lag = 1, k = 1, axis = 2) {
  x <- as.Constant(x)
  if((axis == 2 && ndim(x) < 2) || ndim(x) == 0)
    stop("Invalid axis given input dimensions.")
  else if(axis == 1)
    x <- t(x)
  
  if(ndim(x) == 1)
    m <- size(x)
  else
    m <- dim(x)[setdiff(seq_len(ndim(x)), axis)]
  
  if(k < 0 || k >= m)
    stop("Must have k >= 0 and x must have < k elements along collapsed axis.")
  if(lag <= 0 || lag >= m)
    stop("Must have lag > 0 and x must have < lag elements along collapsed axis.")

  d <- x
  for(i in seq_len(k)) {
    if(ndim(x) == 2)
      d <- d[(1+lag):m,] - d[1:(m-lag),]
    else
      d <- d[(1+lag):m] - d[1:(m-lag)]
    m <- m-1
  }

  if(axis == 1)
    t(d)
  else
    d
}

#'
#' The HStack class.
#'
#' Horizontal concatenation of values.
#'
#' @slot ... \linkS4class{Expression} objects or matrices. All arguments must have the same dimensions except for axis 2 (columns).
#' @name HStack-class
#' @aliases HStack
#' @rdname HStack-class
.HStack <- setClass("HStack", contains = "AffAtom")

#' @param ... \linkS4class{Expression} objects or matrices. All arguments must have the same dimensions except for axis 2 (columns).
#' @rdname HStack-class
HStack <- function(...) {
  arg_list <- lapply(list(...), as.Constant)
  for(idx in seq_along(arg_list)) {
    arg <- arg_list[[idx]]
    if(ndim(arg) == 0)
      arg_list[[idx]] <- flatten(arg)
  }
  .HStack(atom_args = arg_list)
}

#' @param object A \linkS4class{HStack} object.
#' @param values A list of arguments to the atom.
#' @describeIn HStack Horizontally concatenate the values using \code{cbind}.
setMethod("to_numeric", "HStack", function(object, values) {
  # do.call("cbind", values)   # Doesn't work on some objects like xts.
  mat <- Reduce("cbind", values)
  colnames(mat) <- NULL   # Get rid of init column name.
  return(mat)
})

#' @describeIn HStack The dimensions of the atom.
setMethod("dim_from_args", "HStack", function(object) {
  if(ndim(object@args[[1]]) == 1)
    # return(c(sum(sapply(object@args, size)), NA))
    return(c(sum(sapply(object@args, size)), 1))
  else {
    cols <- sum(sapply(object@args, function(arg) { dim(arg)[2] }))
    arg_dim <- dim(object@args[[1]])
    dims <- c(arg_dim[1], cols)
    if(length(arg_dim) >= 3)
      dims <- c(dims, arg_dim[3:length(arg_dim)])
    return(dims)
  }
})

#' @describeIn HStack Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "HStack", function(object) { TRUE })

#' @describeIn HStack Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "HStack", function(object) { TRUE })

#' @describeIn HStack Check all arguments have the same height.
setMethod("validate_args", "HStack", function(object) {
  model <- dim(object@args[[1]])
  error <- "All the input dimensions except for axis 2 (columns) must match exactly."
  len <- length(object@args)
  
  if(len >= 2) {
    for(arg in object@args[2:len]) {
      if(length(dim(arg)) != length(model))
        stop(error)
      else if(length(model) > 1) {
        for(i in 1:length(model)) {
          if(i != 2 && dim(arg)[i] != model[i])
            stop(error)
        }
      }
    }
  }
})

HStack.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.hstack(arg_objs, dim), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn HStack The graph implementation of the atom.
setMethod("graph_implementation", "HStack", function(object, arg_objs, dim, data = NA_real_) {
  HStack.graph_implementation(arg_objs, dim, data)
})

setMethod("cbind2", signature(x = "Expression", y = "ANY"), function(x, y, ...) { HStack(x, y) })
setMethod("cbind2", signature(x = "ANY", y = "Expression"), function(x, y, ...) { HStack(x, y) })

#'
#' The Imag class.
#'
#' This class represents the imaginary part of an expression.
#'
#' @slot expr An \linkS4class{Expression} representing a vector or matrix.
#' @name Imag-class
#' @aliases Imag
#' @rdname Imag-class
.Imag <- setClass("Imag", representation(expr = "Expression"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} representing a vector or matrix.
#' @rdname Imag-class
Imag <- function(expr) { .Imag(expr = expr) }

setMethod("initialize", "Imag", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object An \linkS4class{Imag} object.
#' @param values A list of arguments to the atom.
#' @describeIn Imag The imaginary part of the given value.
setMethod("to_numeric", "Imag", function(object, values) { Im(values[[1]]) })

#' @describeIn Imag The dimensions of the atom.
setMethod("dim_from_args", "Imag", function(object) { dim(object@args[[1]]) })

#' @describeIn Imag Is the atom imaginary?
setMethod("is_imag", "Imag", function(object) { FALSE })

#' @describeIn Imag Is the atom complex valued?
setMethod("is_complex", "Imag", function(object) { FALSE })

#' @describeIn Imag Is the atom symmetric?
setMethod("is_symmetric", "Imag", function(object) { is_hermitian(object@args[[1]]) })

#'
#' The Index class.
#'
#' This class represents indexing or slicing into a matrix.
#'
#' @slot expr An \linkS4class{Expression} representing a vector or matrix.
#' @slot key A list containing the start index, end index, and step size of the slice.
#' @name Index-class
#' @aliases Index
#' @rdname Index-class
.Index <- setClass("Index", representation(expr = "Expression", key = "list"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} representing a vector or matrix.
#' @param key A list containing the start index, end index, and step size of the slice.
#' @rdname Index-class
Index <- function(expr, key) { .Index(expr = expr, key = key) }

setMethod("initialize", "Index", function(.Object, ..., expr, key) {
  .Object@key <- ku_validate_key(key, dim(expr))   # TODO: Double check key validation
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object An \linkS4class{Index} object.
#' @param values A list of arguments to the atom.
#' @describeIn Index The index/slice into the given value.
setMethod("to_numeric", "Index", function(object, values) {
  ku_slice_mat(values[[1]], object@key)
})

#' @describeIn Index The dimensions of the atom.
setMethod("dim_from_args", "Index", function(object) {
  ku_dim(object@key, dim(object@args[[1]]))
})

#' @describeIn Index Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Index", function(object) { TRUE })

#' @describeIn Index Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Index", function(object) { TRUE })

#' @describeIn Index A list containing \code{key}.
setMethod("get_data", "Index", function(object) { list(object@key) })

Index.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  obj <- lo.index(arg_objs[[1]], dim, data[[1]])
  list(obj, list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Index The graph implementation of the atom.
setMethod("graph_implementation", "Index", function(object, arg_objs, dim, data = NA_real_) {
  Index.graph_implementation(arg_objs, dim, data)
})

#'
#' The SpecialIndex class.
#'
#' This class represents indexing using logical indexing or a list of indices into a matrix.
#'
#' @slot expr An \linkS4class{Expression} representing a vector or matrix.
#' @slot key A list containing the start index, end index, and step size of the slice.
#' @name SpecialIndex-class
#' @aliases SpecialIndex
#' @rdname SpecialIndex-class
.SpecialIndex <- setClass("SpecialIndex", representation(expr = "Expression", key = "list", .select_mat = "ConstVal", .dim = "NumORNULL"),
                                          prototype(.select_mat = NA_real_, .dim = NA_real_), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} representing a vector or matrix.
#' @param key A list containing the start index, end index, and step size of the slice.
#' @rdname SpecialIndex-class
SpecialIndex <- function(expr, key) { .SpecialIndex(expr = expr, key = key) }

setMethod("initialize", "SpecialIndex", function(.Object, ..., expr, key) {
  .Object@expr <- expr
  .Object@key <- key
  row <- key[[1]]
  col <- key[[2]]
  
  # Order the entries of expr and select them using key.
  expr_dim <- dim(expr)
  expr_size <- size(expr)
  
  idx_mat <- seq(expr_size)
  idx_mat <- matrix(idx_mat, nrow = expr_dim[1], ncol = expr_dim[2])
  if(is.matrix(row) && is.null(col))
    select_mat <- matrix(idx_mat[row], ncol = 1)
  else if(is.null(row) && is.null(col))
    select_mat <- idx_mat
  else if(is.null(row) && !is.null(col))
    select_mat <- idx_mat[ , col, drop = FALSE]
  else if(!is.null(row) && is.null(col))
    select_mat <- idx_mat[row, , drop = FALSE]
  else
    select_mat <- idx_mat[row, col, drop = FALSE]
  
  .Object@.select_mat <- select_mat
  .Object@.dim <- dim(.Object@.select_mat)
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param x,object An \linkS4class{Index} object.
#' @describeIn SpecialIndex Returns the index in string form.
setMethod("name", "SpecialIndex", function(x) { paste(name(x@args[[1]]), as.character(x@key)) })

#' @param values A list of arguments to the atom.
#' @describeIn Index The index/slice into the given value.
setMethod("to_numeric", "SpecialIndex", function(object, values) {
  ku_slice_mat(values[[1]], object@key)
})

#' @describeIn Index The dimensions of the atom.
setMethod("dim_from_args", "SpecialIndex", function(object) { object@.dim })

#' @describeIn SpecialIndex Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "SpecialIndex", function(object) { TRUE })

#' @describeIn SpecialIndex Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "SpecialIndex", function(object) { TRUE })

#' @describeIn SpecialIndex A list containing \code{key}.
setMethod("get_data", "SpecialIndex", function(object) { list(object@key) })

#' @describeIn SpecialIndex Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "SpecialIndex", function(object) {
  select_mat <- object@.select_mat
  
  if(!is.null(dim(select_mat)))
    final_dim <- dim(select_mat)
  else   # Always cast 1-D arrays as column vectors
    final_dim <- c(length(select_mat), 1)
  
  # Select the chosen entries from expr.
  select_vec <- as.vector(select_mat)
  ##select_vec <- as.matrix(select_mat, nrow=final_dim[1L], ncol=final_dim[2L])
  
  expr_size <- size(object@args[[1]])
  identity <- sparseMatrix(i = 1:expr_size, j = 1:expr_size, x = rep(1, expr_size))
  idmat <- matrix(identity[select_vec, ], ncol = expr_size)
  v <- Vec(expr)
  if(is_scalar(v) || is_scalar(as.Constant(idmat)))
    lowered <- Reshape(idmat * v, c(final_dim[1], final_dim[2]))
  else
    lowered <- Reshape(idmat %*% v, c(final_dim[1], final_dim[2]))
  return(grad(lowered))
})

SpecialIndex.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  select_mat <- object@.select_mat
  final_dim <- dim(object@.select_mat)
  select_vec <- matrix(select_mat, prod(dim(select_mat)), byrow = FALSE)
  
  # Select the chosen entries from expr.
  arg <- arg_objs[[1]]
  arg_size <- size(object@args[[1]])
  id_mat <- sparseMatrix(i = seq_len(arg_size), j = seq_len(arg_size), x = rep(1, arg_size))
  vec_arg <- lu_reshape(arg, c(arg_size, 1))
  mul_mat <- id_mat[select_vec]
  mul_const <- lu_create_const(mul_mat, dim(mul_mat), sparse = TRUE)
  mul_expr <- lu_mul_expr(mul_const, vec_arg, c(nrow(mul_mat), 1))
  obj <- lu_reshape(mul_expr, final_dim)
  return(list(obj, list()))
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Index The graph implementation of the atom.
setMethod("graph_implementation", "SpecialIndex", function(object, arg_objs, dim, data = NA_real_) {
  SpecialIndex.graph_implementation(arg_objs, dim, data)
})

#'
#' The Kron class.
#'
#' This class represents the kronecker product.
#'
#' @slot lh_exp An \linkS4class{Expression} or numeric constant representing the left-hand matrix.
#' @slot rh_exp An \linkS4class{Expression} or numeric constant representing the right-hand matrix.
#' @name Kron-class
#' @aliases Kron
#' @rdname Kron-class
.Kron <- setClass("Kron", representation(lh_exp = "ConstValORExpr", rh_exp = "ConstValORExpr"), contains = "AffAtom")

#' @param lh_exp An \linkS4class{Expression} or numeric constant representing the left-hand matrix.
#' @param rh_exp An \linkS4class{Expression} or numeric constant representing the right-hand matrix.
#' @rdname Kron-class
Kron <- function(lh_exp, rh_exp) { .Kron(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("initialize", "Kron", function(.Object, ..., lh_exp, rh_exp) {
  .Object@lh_exp <- lh_exp
  .Object@rh_exp <- rh_exp
  callNextMethod(.Object, ..., atom_args = list(.Object@lh_exp, .Object@rh_exp))
})

#' @param object A \linkS4class{Kron} object.
#' @param values A list of arguments to the atom.
#' @describeIn Kron The kronecker product of the two values.
setMethod("to_numeric", "Kron", function(object, values) {
  base::kronecker(values[[1]], values[[2]])
})

#' @describeIn Kron Check both arguments are vectors and the first is a constant.
setMethod("validate_args", "Kron", function(object) {
  if(!(is_constant(object@args[[1]]) || is_constant(object@args[[2]])))
    stop("The first argument to Kron must be constant.")
  else if(ndim(object@args[[1]]) != 2 || ndim(object@args[[2]]) != 2)
    stop("Kron requires matrix arguments.")
})

#' @describeIn Kron The dimensions of the atom.
setMethod("dim_from_args", "Kron", function(object) {
  rows <- dim(object@args[[1]])[1] * dim(object@args[[2]])[1]
  cols <- dim(object@args[[1]])[2] * dim(object@args[[2]])[2]
  return(c(rows, cols))
})

#' @describeIn Kron Is the atom convex?
setMethod("is_atom_convex", "Kron", function(object) {
  if(dpp_scope_active()) {
    # Kron is not DPP if any parameters are present
    x <- object@args[[1]]
    y <- object@args[[2]]
    return((is_constant(x) || is_constant(y)) && (is_param_free(x) && is_param_free(y)))
  } else
    return(is_constant(object@args[[1]]) || is_constant(object@args[[2]]))
})

#' @describeIn Kron Is the atom concave?
setMethod("is_atom_concave", "Kron", function(object) {
  return(is_atom_convex(object))
})

#' @describeIn Kron The sign of the atom.
setMethod("sign_from_args", "Kron", function(object) { mul_sign(object@args[[1]], object@args[[2]]) })

#' @param idx An index into the atom.
#' @describeIn Kron Is the composition non-decreasing in argument \code{idx}?
setMethod("is_incr", "Kron", function(object, idx) {
  cst_loc <- ifelse(is_constant(object@args[[1]]), 1, 2)
  is_nonneg(object@args[[cst_loc]])
})

#' @describeIn Kron Is the composition non-increasing in argument \code{idx}?
setMethod("is_decr", "Kron", function(object, idx) {
  cst_loc <- ifelse(is_constant(object@args[[1]]), 1, 2)
  is_nonpos(object@args[[1]])
})

#' @describeIn Kron Is the atom a positive semidefinite matrix?
setMethod("is_psd", "Kron", function(object) {
  # Check a *sufficient condition* that the expression is PSD, 
  # by checking if both arguments are PSD or both are NSD.
  case1 <- is_psd(object@args[[1]]) && is_psd(object@args[[2]])
  case2 <- is_nsd(object@args[[1]]) && is_nsd(object@args[[2]])
  return(case1 || case2)
})

#' @describeIn Kron Is the atom a negative semidefinite matrix?
setMethod("is_nsd", "Kron", function(object) {
  # Check a *sufficient condition* that the expression is NSD, 
  # by checking if one argument is PSD and the other is NSD.
  case1 <- is_psd(object@args[[1]]) && is_nsd(object@args[[2]])
  case2 <- is_nsd(object@args[[1]]) && is_psd(object@args[[2]])
  return(case1 || case2)
})

Kron.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  if(is_constant(object@args[[1]]))
    list(lo.kron_r(arg_objs[[1]], arg_objs[[2]], dim), list())
  else
    list(lo.kron_l(arg_objs[[1]], arg_objs[[2]], dim), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Kron The graph implementation of the atom.
setMethod("graph_implementation", "Kron", function(object, arg_objs, dim, data = NA_real_) {
  Kron.graph_implementation(arg_objs, dim, data)
})

#'
#' The PartialTrace atom.
#' 
#' Assume \eqn{expr = X_1 \otimes \cdots \otimes X_n} is a 2D Kronecker product 
#' composed of \eqn{n = length(dims)} implicit subsystems. 
#' Letting \eqn{k = axis}, this class represents the partial trace of \eqn{expr} 
#' along its \eqn{k}-th implicit subsystem:
#' 
#' \deqn{tr(X_k) (X_1 \otimes \cdots \otimes X_{k-1} \otimes X_{k+1} \otimes \cdots \otimes X_n)}
#'
#' @param expr An \linkS4class{Expression} representing a 2D expression of which to take the partial trace.
#' @param dims A vector of integers encoding the dimensions of each subsystem.
#' @param axis The index of the subsystem to be traced out from the tensor product that defines \code{expr}.
#' @return The partial trace of \code{expr}.
PartialTrace <- function(expr, dims, axis) {
  expr <- as.Constant(expr)
  dims <- as.integer(dims)
  axis <- as.integer(axis)
  
  if(any(dims) <= 0)
    stop("dims must have positive integer entries")
  if(axis <= 0)
    stop("axis must be a positive integer")
  if(ndim(expr) < 2 || dim(expr)[1] != dim(expr)[2])
    stop("Only supports square matrices")
  if(axis <= 0 || axis > length(dims))
    stop("Invalid axis argument, should be between 1 and ", length(dims))
  if(dim(expr)[1] != base::prod(dims))
    stop("Dimension of system doesn't correspond to dimension of subsystems")
  
  term_list <- lapply(seq_len(dims[axis]), function(j) { PartialTrace.term(expr, j, dims, axis) })
  return(AddExpression(arg_groups = term_list))
}

#
# Helper function for PartialTrace.
#
# Parameters
# -----------
# expr: The 2D expression of which to take the partial trace.
# j: The term in the partial trace sum.
# dims: A vector of integers encoding the dimensions of each subsystem.
# axis: The index of the subsystem to be traced out from the tensor product that defines expr.
#
# (I ⊗ <j| ⊗ I) x (I ⊗ |j> ⊗ I) for all j's
# in the system we want to trace out.
# This function returns the jth term in the sum, namely
# (I ⊗ <j| ⊗ I) x (I ⊗ |j> ⊗ I).
#
PartialTrace.term <- function(expr, j, dims, axis) {
  a <- Matrix(1, 1, 1, sparse = TRUE)
  b <- Matrix(1, 1, 1, sparse = TRUE)
  
  ndims <- length(dims)
  for(i_axis in seq_len(ndims)) {
    dim <- dims[i_axis]
    if(i_axis == axis) {
      v <- sparseMatrix(j, 1, x = 1, dims = c(dim, 1))
      a <- kronecker(a, t(v))
      b <- kronecker(b, v)
    } else {
      eye_mat <- sparseMatrix(seq_len(dim), seq_len(dim), x = 1)
      a <- kronecker(a, eye_mat)
      b <- kronecker(b, eye_mat)
    }
  }
  return(a %*% expr %*% b)
}

#'
#' The PartialTranspose atom.
#' 
#' Assume \eqn{expr = X_1 \otimes \cdots \otimes X_n} is a 2D Kronecker product 
#' composed of \eqn{n = length(dims)} implicit subsystems. 
#' Letting \eqn{k = axis}, this class represents the partial transpose of \eqn{expr}, 
#' with the transpose applied to its \eqn{k}-th implicit subsystem:
#' 
#' \deqn{X_1 \otimes ...\otimes X_k^T \otimes ... \otimes X_n}.
#'
#' @param expr An \linkS4class{Expression} representing a 2D expression of which to take the partial transpose.
#' @param dims A vector of integers encoding the dimensions of each subsystem.
#' @param axis The index of the subsystem to be transposed out from the tensor product that defines \code{expr}.
#' @return The partial trace of \code{expr}.
PartialTrace <- function(expr, dims, axis) {
  expr <- as.Constant(expr)
  dims <- as.integer(dims)
  axis <- as.integer(axis)
  
  if(any(dims) <= 0)
    stop("dims must have positive integer entries")
  if(axis <= 0)
    stop("axis must be a positive integer")
  if(ndim(expr) < 2 || dim(expr)[1] != dim(expr)[2])
    stop("Only supports square matrices")
  if(axis <= 0 || axis > length(dims))
    stop("Invalid axis argument, should be between 1 and ", length(dims))
  if(dim(expr)[1] != base::prod(dims))
    stop("Dimension of system doesn't correspond to dimension of subsystems")
  
  term_list <- list()
  for(i in seq_len(dims[axis])) {
    for(j in seq_len(dims[axis]))
      term_list <- c(term_list, PartialTranspose.term(expr, i, j, dims, axis))
  }
  return(AddExpression(arg_groups = term_list))
}

#
# Helper function for PartialTranspose.
#
# Parameters
# -----------
# expr: The 2D expression of which to take the partial transpose.
# i: The term in the partial transpose sum.
# j: The term in the partial transpose sum.
# dims: A vector of integers encoding the dimensions of each subsystem.
# axis: The index of the subsystem to be transposed from the tensor product that defines expr.
#
# (I ⊗ |i><j| ⊗ I) x (I ⊗ |i><j| ⊗ I) for all (i,j)'s
# in the system we want to transpose.
# This function returns the (i,j)-th term in the sum, namely
# (I ⊗ |i><j| ⊗ I) x (I ⊗ |i><j| ⊗ I).
#
PartialTranspose.term <- function(expr, i, j, dims, axis) {
  a <- Matrix(1, 1, 1, sparse = TRUE)
  
  ndims <- length(dims)
  for(i_axis in seq_len(ndims)) {
    dim <- dims[i_axis]
    if(i_axis == axis) {
      v <- sparseMatrix(i, j, x = 1, dims = c(dim, dim))
      a <- kronecker(a, v)
    } else {
      eye_mat <- sparseMatrix(seq_len(dim), seq_len(dim), x = 1)
      a <- kronecker(a, eye_mat)
    }
  }
  return(a %*% expr %*% a)
}

#'
#' The Promote class.
#' 
#' This class represents the promotion of a scalar expression into a vector/matrix.
#'
#' @slot expr An \linkS4class{Expression} or numeric constant.
#' @slot promoted_dim The desired dimensions.
#' @name Promote-class
#' @aliases Promote
#' @rdname Promote-class
.Promote <- setClass("Promote", representation(expr = "Expression", promoted_dim = "numeric"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} or numeric constant.
#' @param promoted_dim The desired dimensions.
#' @rdname Promote-class
Promote <- function(expr, promoted_dim) { .Promote(expr = expr, promoted_dim = promoted_dim) } 

promote <- function(expr, promoted_dim) { 
  expr <- as.Constant(expr)
  if(!all(dim(expr) == promoted_dim)) {
    if(!is_scalar(expr))
      stop("Only scalars may be promoted.")
    return(Promote(expr = expr, promoted_dim = promoted_dim))
  } else
    return(expr)
}

setMethod("initialize", "Promote", function(.Object, ..., expr, promoted_dim) {
  .Object@expr <- expr
  .Object@promoted_dim <- promoted_dim
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{Promote} object.
#' @param values A list containing the value to promote.
#' @describeIn Promote Promotes the value to the new dimensions.
setMethod("to_numeric", "Promote", function(object, values) {
  array(1, dim = object@promoted_dim) * as.vector(values[[1]])[1]
})

#' @describeIn Promote Is the expression symmetric?
setMethod("is_symmetric", "Promote", function(object) {
  ndim(object) == 2 && dim(object)[1] == dim(object)[2]
})

#' @describeIn Promote Returns the (row, col) dimensions of the expression.
setMethod("dim_from_args", "Promote", function(object) { object@promoted_dim })

#' @describeIn Promote Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Promote", function(object) { TRUE })

#' @describeIn Promote Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Promote", function(object) { TRUE })

#' @describeIn Promote Returns information needed to reconstruct the expression besides the args.
setMethod("get_data", "Promote", function(object) { list(object@promoted_dim) })

Promote.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.promote(arg_objs[[1]], dim), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Promote The graph implementation of the atom.
setMethod("graph_implementation", "Promote", function(object, arg_objs, dim, data = NA_real_) {
  Promote.graph_implementation(arg_objs, dim, data)
})

#'
#' The Real class.
#'
#' This class represents the real part of an expression.
#'
#' @slot expr An \linkS4class{Expression} representing a vector or matrix.
#' @name Real-class
#' @aliases Real
#' @rdname Real-class
.Real <- setClass("Real", representation(expr = "Expression"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} representing a vector or matrix.
#' @rdname Real-class
Real <- function(expr) { .Real(expr = expr) }

setMethod("initialize", "Real", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object An \linkS4class{Real} object.
#' @param values A list of arguments to the atom.
#' @describeIn Real The imaginary part of the given value.
setMethod("to_numeric", "Real", function(object, values) { Re(values[[1]]) })

#' @describeIn Real The dimensions of the atom.
setMethod("dim_from_args", "Real", function(object) { dim(object@args[[1]]) })

#' @describeIn Real Is the atom imaginary?
setMethod("is_imag", "Real", function(object) { FALSE })

#' @describeIn Real Is the atom complex valued?
setMethod("is_complex", "Real", function(object) { FALSE })

#' @describeIn Real Is the atom symmetric?
setMethod("is_symmetric", "Real", function(object) { is_hermitian(object@args[[1]]) })

#'
#' The Reshape class.
#'
#' This class represents the reshaping of an expression. The operator vectorizes the expression,
#' then unvectorizes it into the new dimensions. Entries are stored in column-major order.
#'
#' @slot expr An \linkS4class{Expression} or numeric matrix.
#' @slot new_dim The new dimensions.
#' @name Reshape-class
#' @aliases Reshape
#' @rdname Reshape-class
.Reshape <- setClass("Reshape", representation(expr = "ConstValORExpr", new_dim = "numeric", byrow = "logical"), 
                                prototype(byrow = FALSE), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} or numeric matrix.
#' @param new_dim The new dimensions.
#' @param byrow A logical value indicating whether the matrix is filled by rows. If \code{FALSE} (default), the matrix is filled by columns.
#' @rdname Reshape-class
Reshape <- function(expr, new_dim, byrow = FALSE) { .Reshape(expr = expr, new_dim = new_dim, byrow = byrow) }

setMethod("initialize", "Reshape", function(.Object, ..., expr, new_dim, byrow) {
  if(length(new_dim) > 2)
    stop("Expressions of dimension greater than 2 are not supported")
  .Object@new_dim <- new_dim
  .Object@expr <- expr
  .Object@byrow <- byrow
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{Reshape} object.
#' @param values A list of arguments to the atom.
#' @describeIn Reshape Reshape the value into the specified dimensions.
setMethod("to_numeric", "Reshape", function(object, values) {
  if(object@byrow) {
    mat <- values[[1]]
    dim(mat) <- rev(object@new_dim)
    return(t(mat))
  } else {
    dim(values[[1]]) <- object@new_dim
    return(values[[1]])
  }

})

#' @describeIn Reshape Check the new shape has the same number of entries as the old.
setMethod("validate_args", "Reshape", function(object) {
  old_len <- size(object@args[[1]])
  new_len <- prod(object@new_dim)
  if(old_len != new_len)
    stop("Invalid reshape dimensions (", paste(object@new_dim, sep = ","), ")")
})

#' @describeIn Reshape The \code{c(rows, cols)} dimensions of the new expression.
setMethod("dim_from_args", "Reshape", function(object) { object@new_dim })

#' @describeIn Reshape Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Reshape", function(object) { TRUE })

#' @describeIn Reshape Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Reshape", function(object) { TRUE })

#' @describeIn Reshape Returns a list containing the new shape.
setMethod("get_data", "Reshape", function(object) { list(object@new_dim, object@byrow) })

Reshape.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  arg <- arg_objs[[1]]
  if(data[2]) {   # byrow = TRUE
    arg <- lo.transpose(arg)
    if(length(dim) <= 1)
      list(lo.reshape(arg, dim), list())
    else {
      result <- lo.reshape(arg, c(dim[2], dim[1]))
      list(lo.transpose(result), list())
    }
  } else   # byrow = FALSE
    list(lo.reshape(arg, dim), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Reshape The graph implementation of the atom.
setMethod("graph_implementation", "Reshape", function(object, arg_objs, dim, data = NA_real_) {
  Reshape.graph_implementation(arg_objs, dim, data)
})

Reshape.deep_flatten <- function(x)
{
  # Base cases
  if(is(x, "Expression")) {
    if(length(dim(x)) == 1 || all(dim(x) == c(1,1)))
      return(x)
    else
      return(flatten(x))
  } else if(is.numeric(x)) {
    x <- as.Constant(x)
    return(flatten(x))
  }
  
  # Recursion
  if(is(x, "list")) {
    y <- list()
    for(xo in x) {
      x1 <- deep_flatten(x0)
      y <- c(y, x1)
    }
    y <- HStack(y)
    return(y)
  }
  
  stop("The input to deep_flatten must be an Expression, numeric array or scalar, or a nested list thereof. Received input of type ", class(x))
}

#'
#' The SumEntries class.
#'
#' This class represents the sum of all entries in a vector or matrix.
#'
#' @slot expr An \linkS4class{Expression} representing a vector or matrix.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @slot keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @name SumEntries-class
#' @aliases SumEntries
#' @rdname SumEntries-class
.SumEntries <- setClass("SumEntries", contains = c("AxisAtom", "AffAtom"))

#' @param expr An \linkS4class{Expression} representing a vector or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @rdname SumEntries-class
SumEntries <- function(expr, axis = NA_real_, keepdims = FALSE) { .SumEntries(expr = expr, axis = axis, keepdims = keepdims) }

#' @param object A \linkS4class{SumEntries} object.
#' @param values A list of arguments to the atom.
#' @describeIn SumEntries Sum the entries along the specified axis.
setMethod("to_numeric", "SumEntries", function(object, values) {
  apply_with_keepdims(values[[1]], sum, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn SumEntries Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "SumEntries", function(object) { TRUE })

#' @describeIn SumEntries Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "SumEntries", function(object) { FALSE })

SumEntries.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  # TODO: Handle keepdims properly by setting a 1-D vector's dimension to c(len, NA_integer_).
  axis <- data[[1]]
  keepdims <- data[[2]]
  if(is.na(axis))
    obj <- lo.sum_entries(arg_objs[[1]], dim)
  else if(axis == 1) {
    # if(keepdims)
    #  const_dim <- c(arg_objs[[1]]$dim[2], 1)
    # else
    #  const_dim <- c(arg_objs[[1]]$dim[2], NA_integer_)
    
    # Always treat result as a column vector.
    const_dim <- c(arg_objs[[1]]$dim[2], 1)
    ones <- create_const(array(1, dim = const_dim), const_dim)
    obj <- lo.rmul_expr(arg_objs[[1]], ones, dim)
  } else {   # axis == 2
    # if(keepdims)
    #  const_dim <- c(1, arg_objs[[1]]$dim[1])
    # else
    #  const_dim <- c(arg_objs[[1]]$dim[1], NA_integer_)
    # ones <- create_const(array(1, dim = const_dim), const_dim)
    # obj <- lo.mul_expr(ones, arg_objs[[1]], dim)
    
    if(keepdims) {
      # Keep result as a row vector.
      const_dim <- c(1, arg_objs[[1]]$dim[1])
      ones <- create_const(array(1, dim = const_dim), const_dim)
      obj <- lo.mul_expr(ones, arg_objs[[1]], dim)
    } else {
      # Treat collapsed 1-D vector as a column vector.
      const_dim <- c(arg_objs[[1]]$dim[1], 1)
      ones <- create_const(array(1, dim = const_dim), const_dim)
      obj <- lo.rmul_expr(lo.transpose(arg_objs[[1]]), ones, dim)
    }
  }
  list(obj, list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn SumEntries The graph implementation of the atom.
setMethod("graph_implementation", "SumEntries", function(object, arg_objs, dim, data = NA_real_) {
  SumEntries.graph_implementation(arg_objs, dim, data)
})

#'
#' The Trace class.
#'
#' This class represents the sum of the diagonal entries in a matrix.
#'
#' @slot expr An \linkS4class{Expression} representing a matrix.
#' @name Trace-class
#' @aliases Trace
#' @rdname Trace-class
.Trace <- setClass("Trace", representation(expr = "Expression"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} representing a matrix.
#' @rdname Trace-class
Trace <- function(expr) { .Trace(expr = expr) }

setMethod("initialize", "Trace", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{Trace} object.
#' @param values A list of arguments to the atom.
#' @describeIn Trace Sum the diagonal entries.
setMethod("to_numeric", "Trace", function(object, values) { sum(diag(values[[1]])) })

#' @describeIn Trace Check the argument is a square matrix.
setMethod("validate_args", "Trace", function(object) {
  arg_dim <- dim(object@args[[1]])
  if(arg_dim[1] != arg_dim[2])
    stop("Argument to Trace must be a square matrix")
})

#' @describeIn Trace The atom is a scalar.
setMethod("dim_from_args", "Trace", function(object){ c(1,1) })

#' @describeIn Trace The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "Trace", function(object) {
  is_nonneg <- is_nonneg(object@args[[1]]) || is_psd(object@args[[1]])
  is_nonpos <- is_nonpos(object@args[[1]]) || is_nsd(object@args[[1]])
  c(is_nonneg, is_nonpos)
})

#' @describeIn Trace A logical value indicating whether the atom is real.
setMethod("is_real", "Trace", function(object) {
  is_real(object@args[[1]]) || is_hermitian(object@args[[1]])
})

#' @describeIn Trace A logical value indicating whether the atom is complex.
setMethod("is_complex", "Trace", function(object) { !is_real(object) })

#' @describeIn Trace Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Trace", function(object) { TRUE })

#' @describeIn Trace Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Trace", function(object) { FALSE })

Trace.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.trace(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Trace The graph implementation of the atom.
setMethod("graph_implementation", "Trace", function(object, arg_objs, dim, data = NA_real_) {
  Trace.graph_implementation(arg_objs, dim, data)
})

#'
#' The Transpose class.
#'
#' This class represents the matrix transpose.
#'
#' @name Transpose-class
#' @aliases Transpose
#' @rdname Transpose-class
.Transpose <- setClass("Transpose", representation(expr = "Expression", axes = "ConstValORNULL"), prototype(axes = NULL), contains = "AffAtom")

Transpose <- function(expr, axes = NULL) { .Transpose(expr = expr, axes = axes) }

setMethod("initialize", "Transpose", function(.Object, ..., expr, axes = NULL) {
  .Object@expr <- expr
  .Object@axes <- axes
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{Transpose} object.
#' @param values A list of arguments to the atom.
#' @describeIn Transpose The transpose of the given value.
setMethod("to_numeric", "Transpose", function(object, values) {
  if(is.vector(values[[1]]))
    return(t(values[[1]]))
  else if(is(values[[1]], "Matrix")) {
    if(!is.null(object@axes))
      stop("Cannot permute Matrix object axes to (", paste(object@axes, collapse = ","), ")")
    return(t(values[[1]]))
  } else
    return(aperm(values[[1]], perm = object@axes))
})

#' @describeIn Transpose Is the expression symmetric?
setMethod("is_symmetric", "Transpose", function(object) { is_symmetric(object@args[[1]]) })

#' @describeIn Transpose Is the expression skew symmetric?
setMethod("is_skew_symmetric", "Transpose", function(object) { is_skew_symmetric(object@args[[1]]) })

#' @describeIn Transpose Is the expression hermitian?
setMethod("is_hermitian", "Transpose", function(object) { is_hermitian(object@args[[1]]) })

#' @describeIn Transpose The dimensions of the atom.
setMethod("dim_from_args", "Transpose", function(object) {
  if(is.null(object@axes))
    rev(dim(object@args[[1]]))
  else
    dim(object@args[[1]])[object@axes]
})

#' @describeIn Transpose Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Transpose", function(object) { TRUE })

#' @describeIn Transpose Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Transpose", function(object) { TRUE })

#' @describeIn Transpose Returns the axes for transposition.
setMethod("get_data", "Transpose", function(object) { list(object@axes) })

Transpose.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.transpose(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Transpose The graph implementation of the atom.
setMethod("graph_implementation", "Transpose", function(object, arg_objs, dim, data = NA_real_) {
  Transpose.graph_implementation(arg_objs, dim, data)
})

#'
#' The UpperTri class.
#'
#' The vectorized strictly upper triagonal entries of a matrix.
#'
#' @slot expr An \linkS4class{Expression} or numeric matrix.
#' @name UpperTri-class
#' @aliases UpperTri
#' @rdname UpperTri-class
.UpperTri <- setClass("UpperTri", representation(expr = "ConstValORExpr"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} or numeric matrix.
#' @rdname UpperTri-class
UpperTri <- function(expr) { .UpperTri(expr = expr) }

setMethod("initialize", "UpperTri", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object An \linkS4class{UpperTri} object.
#' @param values A list of arguments to the atom.
#' @describeIn UpperTri Vectorize the upper triagonal entries.
setMethod("to_numeric", "UpperTri", function(object, values) {
  # Vectorize the upper triagonal entries
  tridx <- upper.tri(values[[1]], diag = FALSE)
  values[[1]][tridx]
})

#' @describeIn UpperTri Check the argument is a square matrix.
setMethod("validate_args", "UpperTri", function(object) {
  arg_dim <- dim(object@args[[1]])
  if(ndim(object@args[[1]]) != 2 || arg_dim[1] != arg_dim[2])
    stop("Argument to UpperTri must be a square matrix.")
})

#' @describeIn UpperTri The dimensions of the atom.
setMethod("dim_from_args", "UpperTri", function(object) {
  arg_dim <- dim(object@args[[1]])
  rows <- arg_dim[1]
  cols <- arg_dim[2]
  c(floor(rows*(cols-1)/2), 1)
})

#' @describeIn UpperTri Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "UpperTri", function(object) { TRUE })

#' @describeIn UpperTri Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "UpperTri", function(object) { TRUE })

UpperTri.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.upper_tri(arg_objs[[1]]), list())
}

UpperTri.vec_to_upper_tri <- function(expr, strict = FALSE) {
  expr <- as.Constant(expr)
  ell <- dim(expr)[[1]]
  if(strict) {
    # n * (n-1)/2 == ell
    n <- floor(((8*ell + 1)^0.5 + 1)/2)
  } else {
    # n * (n+1)/2 == ell
    n <- floor(((8*ell + 1)^0.5 - 1)/2)
  }
  
  # Form a matrix P of dimensions (n^2, ell).
  #     the i-th block of n rows of P gives the entries of the i-th row
  #     of the upper-triangular matrix associated with expr.
  # Compute expr2 <- P %*% expr
  # Compute expr3 <- t(reshape(expr2, dim = c(n,n)))
  #     expr3 is the matrix formed by reading length-n blocks of expr 2,
  #     and letting each block for a row of expr3.
  P_rows <- c()
  P_row <- 1
  for(mat_row in seq(n)) {
    entries_in_row <- n - mat_row + 1
    if(strict)
      entries_in_row <- entries_in_row - 1
    P_row <- P_row + n - entries_in_row   # these are zeros
    P_rows <- c(P_rows, seq(P_row, P_row + entries_in_row))
    P_row <- P_row + entries_in_row
  }
  P_cols <- seq(ell)
  P_vals <- rep(1, ell)
  P <- sparseMatrix(P_rows, P_cols, x = P_vals, dims = c(n^2, ell))
  expr2 <- P @ expr
  expr3 <- t(Reshape(expr2, c(n, n)))
  return(expr3)
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn UpperTri The graph implementation of the atom.
setMethod("graph_implementation", "UpperTri", function(object, arg_objs, dim, data = NA_real_) {
  UpperTri.graph_implementation(arg_objs, dim, data)
})

# Reshape into single column vector.
Vec <- function(X, byrow = FALSE) {
  X <- as.Constant(X)
  # Reshape(expr = X, new_dim = size(X), byrow = byrow)
  Reshape(expr = X, new_dim = c(size(X), 1), byrow = byrow)
}

#'
#' The VStack class.
#'
#' Vertical concatenation of values.
#'
#' @slot ... \linkS4class{Expression} objects or matrices. All arguments must have the same number of columns.
#' @name VStack-class
#' @aliases VStack
#' @rdname VStack-class
.VStack <- setClass("VStack", contains = "AffAtom")

#' @param ... \linkS4class{Expression} objects or matrices. All arguments must have the same number of columns.
#' @rdname VStack-class
VStack <- function(...) { .VStack(atom_args = list(...)) }

#' @param object A \linkS4class{VStack} object.
#' @param values A list of arguments to the atom.
#' @describeIn VStack Vertically concatenate the values using \code{rbind}.
setMethod("to_numeric", "VStack", function(object, values) {
  # do.call("rbind", values)   # Doesn't work on some objects like xts.
  mat <- Reduce("rbind", values)
  rownames(mat) <- NULL   # Get rid of init row name
  return(mat)
})

#' @describeIn VStack Check all arguments have the same width.
setMethod("validate_args", "VStack", function(object) {
  model <- dim(object@args[[1]])
  if(length(object@args) >= 2) {
    for(arg in object@args[2:length(object@args)]) {
      arg_dim <- dim(arg)
      if(length(arg_dim) != length(model) || length(model) != length(arg_dim) ||
         (length(model) > 1 && any(model[2:length(model)] != arg_dim[2:length(arg_dim)])) ||
         (length(model) <= 1 && any(model != arg_dim)))
        stop("All the input dimensions except for axis 1 must match exactly.")
    }
  }
})

#' @describeIn VStack The dimensions of the atom.
setMethod("dim_from_args", "VStack", function(object) {
  if(ndim(object@args[[1]]) == 0)
    c(length(object@args), 1)
  else if(ndim(object@args[[1]]) == 1)
    c(length(object@args), dim(object@args[[1]])[1])
  else {
    rows <- sum(sapply(object@args, function(arg) { dim(arg)[1] }))
    arg_dim <- dim(object@args[[1]])
    if(length(arg_dim) < 2)
      # c(rows, NA)
      c(rows, 1)
    else
      c(rows, arg_dim[2:length(arg_dim)])
  }
})

#' @describeIn VStack Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "VStack", function(object) { TRUE })

#' @describeIn VStack Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "VStack", function(object) { TRUE })

VStack.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.vstack(arg_objs, dim), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn VStack The graph implementation of the atom.
setMethod("graph_implementation", "VStack", function(object, arg_objs, dim, data = NA_real_) {
  VStack.graph_implementation(arg_objs, dim, data)
})

#'
#' The Wrap class.
#'
#' This virtual class represents a no-op wrapper to assert properties.
#'
#' @name Wrap-class
#' @aliases Wrap
#' @rdname Wrap-class
Wrap <- setClass("Wrap", contains = c("VIRTUAL", "AffAtom"))

#' @param object A \linkS4class{Wrap} object.
#' @param values A list of arguments to the atom.
#' @describeIn Wrap Returns the input value.
setMethod("to_numeric", "Wrap", function(object, values) { values[[1]] })

#' @describeIn Wrap The dimensions of the atom.
setMethod("dim_from_args", "Wrap", function(object) { dim(object@args[[1]]) })

#' @describeIn Wrap Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Wrap", function(object) { TRUE })

#' @describeIn Wrap Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Wrap", function(object) { TRUE })

Wrap.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(arg_objs[[1]], list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Wrap The graph implementation of the atom.
setMethod("graph_implementation", "Wrap", function(object, arg_objs, dim, data = NA_real_) {
  Wrap.graph_implementation(arg_objs, dim, data)
})

#'
#' The PSDWrap class.
#'
#' A no-op wrapper to assert the input argument is positive semidefinite.
#'
#' @name PSDWrap-class
#' @aliases PSDWrap
#' @rdname PSDWrap-class
.PSDWrap <- setClass("PSDWrap", contains = "Wrap")

#' @param arg A \linkS4class{Expression} object or matrix.
#' @rdname PSDWrap-class
PSDWrap <- function(arg) { .PSDWrap(atom_args = list(arg)) }

#' @param object A \linkS4class{PSDWrap} object.
#' @describeIn PSDWrap Check the input is a square matrix.
setMethod("validate_args", "PSDWrap", function(object) {
  arg <- object@args[[1]]
  arg_dim <- dim(arg)
  ndim_test <- (length(arg_dim) == 2)
  if(!ndim_test || arg_dim[1] != arg_dim[2])
    stop("The input must be a square matrix")
})

#' @describeIn PSDWrap Is the atom positive semidefinite?
setMethod("is_psd", "PSDWrap", function(object) { TRUE })

#' @describeIn PSDWrap Is the atom negative semidefinite?
setMethod("is_nsd", "PSDWrap", function(object) { FALSE })

#' @describeIn PSDWrap Is the atom Hermitian?
setMethod("is_hermitian", "PSDWrap", function(object) { TRUE })

validate_real_square <- function(arg) {
  arg_dim <- dim(arg)
  ndim_test <- (length(arg_dim) == 2)
  if(!ndim_test || arg_dim[1] != arg_dim[2])
    stop("The input must be a square matrix")
  else if(!is_real(arg))
    stop("The input must be a real matrix")
}

#'
#' The SymmetricWrap class.
#'
#' A no-op wrapper to assert the input argument is symmetric.
#'
#' @name SymmetricWrap-class
#' @aliases SymmetricWrap
#' @rdname SymmetricWrap-class
.SymmetricWrap <- setClass("SymmetricWrap", contains = "Wrap")

#' @param arg A \linkS4class{Expression} object or matrix.
#' @rdname SymmetricWrap-class
SymmetricWrap <- function(arg) { .SymmetricWrap(atom_args = list(arg)) }

#' @param object A \linkS4class{SymmetricWrap} object.
#' @describeIn SymmetricWrap Check the input is a real square matrix.
setMethod("validate_args", "SymmetricWrap", function(object) {
  validate_real_square(object@args[[1]])
})

#' @describeIn SymmetricWrap Is the atom symmetric?
setMethod("is_symmetric", "SymmetricWrap", function(object) { TRUE })

#' @describeIn SymmetricWrap Is the atom Hermitian?
setMethod("is_hermitian", "SymmetricWrap", function(object) { TRUE })

#'
#' The HermitianWrap class.
#'
#' A no-op wrapper to assert the input argument is Hermitian.
#'
#' @name HermitianWrap-class
#' @aliases HermitianWrap
#' @rdname HermitianWrap-class
.HermitianWrap <- setClass("HermitianWrap", contains = "Wrap")

#' @param arg A \linkS4class{Expression} object or matrix.
#' @rdname HermitianWrap-class
HermitianWrap <- function(arg) { .HermitianWrap(atom_args = list(arg)) }

#' @param object A \linkS4class{HermitianWrap} object.
#' @describeIn HermitianWrap Check the input is a real square matrix.
setMethod("validate_args", "HermitianWrap", function(object) {
  arg <- object@args[[1]]
  arg_dim <- dim(arg)
  ndim_test <- (length(arg_dim) == 2)
  if(!ndim_test || arg_dim[1] != arg_dim[2])
    stop("The input must be a square matrix")
})

#' @describeIn HermitianWrap Is the atom Hermitian?
setMethod("is_hermitian", "HermitianWrap", function(object) { TRUE })

#'
#' The SkewSymmetricWrap class.
#'
#' A no-op wrapper to assert the input argument is skew symmetric.
#'
#' @name SkewSymmetricWrap-class
#' @aliases SkewSymmetricWrap
#' @rdname SkewSymmetricWrap-class
.SkewSymmetricWrap <- setClass("SkewSymmetricWrap", contains = "Wrap")

#' @param arg A \linkS4class{Expression} object or matrix.
#' @rdname SkewSymmetricWrap-class
SkewSymmetricWrap <- function(arg) { .SkewSymmetricWrap(atom_args = list(arg)) }

#' @param object A \linkS4class{SkewSymmetricWrap} object.
#' @describeIn SkewSymmetricWrap Check the input is a real square matrix.
setMethod("validate_args", "SkewSymmetricWrap", function(object) {
  validate_real_square(object@args[[1]])
})

#' @describeIn SkewSymmetricWrap Is the atom skew symmetric?
setMethod("is_skew_symmetric", "SkewSymmetricWrap", function(object) { TRUE })

setMethod("rbind2", signature(x = "Expression", y = "ANY"), function(x, y, ...) { VStack(x, y) })
setMethod("rbind2", signature(x = "ANY", y = "Expression"), function(x, y, ...) { VStack(x, y) })

Bmat <- function(block_lists) {
  row_blocks <- lapply(block_lists, function(blocks) { do.call("HStack", blocks) })
  do.call("VStack", row_blocks)
}

#' 
#' The ScalarProduct atom.
#' 
#' Standard inner product or scalar product of (x,y)
#' 
#' @param x An \linkS4class{Expression} or numeric value. The conjugate-linear argument to the inner product.
#' @param y An \linkS4class{Expression} or numeric value. The linear argument to the inner product.
#' @return An \linkS4class{Expression} representing the standard inner product of (\code{x}, \code{y}), conjugate-linear in \code{x}.
#' @rdname ScalarProduct-int
ScalarProduct <- function(x, y) {
  x <- DeepFlatten(x)
  y <- DeepFlatten(y)
  xy_prod <- Multiply(Conjugate(x), y)
  return(SumEntries(xy_prod))
}
