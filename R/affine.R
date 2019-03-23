#'
#' The AffAtom class.
#'
#' This virtual class represents an affine atomic expression.
#'
#' @name AffAtom-class
#' @aliases AffAtom
#' @rdname AffAtom-class
AffAtom <- setClass("AffAtom", contains = c("VIRTUAL", "Atom"))

#' @describeIn AffAtom Does the atom handle complex numbers?
setMethod("allow_complex", "AffAtom", function(object) { TRUE })

#' @describeIn AffAtom The sign of the atom.
setMethod("sign_from_args", "AffAtom", function(object) { sum_signs(object@args) })

#' @describeIn AffAtom Is the atom imaginary?
setMethod("is_imag", "AffAtom", function(object) { all(sapply(object@args, is_imag)) })

#' @describeIn AffAtom Is the atom complex valued?
setMethod("is_complex", "AffAtom", function(object) { all(sapply(object@args, is_complex)) })

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
  names(var_offsets) <- var_names
  graph <- graph_implementation(object, fake_args, dim(object), get_data(object))
  fake_expr <- graph[[1]]

  # Get the matrix representation of the function.
  prob_mat <- get_problem_matrix(list(fake_expr), var_offsets, NA)
  V <- prob_mat[[1]]
  I <- prob_mat[[2]] + 1   # TODO: R uses 1-indexing, but get_problem_matrix returns with 0-indexing
  J <- prob_mat[[3]] + 1
  dims <- c(offset, size(object))
  stacked_grad <- sparseMatrix(i = J, j = I, x = V, dims = dims)

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

setMethod("name", "AddExpression", function(x) {
  paste(sapply(x@args, as.character), collapse = " + ")
})

#' @param values A list of arguments to the atom.
#' @describeIn AddExpression Sum all the values.
setMethod("to_numeric", "AddExpression", function(object, values) {
  # values <- lapply(values, intf_convert_if_scalar)
  Reduce("+", values)
})

setMethod("is_symmetric", "AddExpression", function(object) {
  symm_args <- all(sapply(object@args, is_symmetric))
  return(dim(object)[1] == dim(object)[2] && symm_args)
})

setMethod("is_hermitian", "AddExpression", function(object) {
  herm_args <- all(sapply(object@args, is_hermitian))
  return(dim(object)[1] == dim(object)[2] && herm_args)
})

# As initialize takes in the arg_groups instead of args, we need a special copy function.
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
UnaryOperator <- setClass("UnaryOperator", representation(expr = "Expression", op_name = "character", op_func = "function"), contains = "AffAtom")

setMethod("initialize", "UnaryOperator", function(.Object, ..., expr, op_name, op_func) {
  .Object@expr <- expr
  .Object@op_name <- op_name
  .Object@op_func <- op_func
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

setMethod("name", "UnaryOperator", function(x) {
  paste(x@op_name, name(x@args[[1]]), sep = "")
})

#' @param object A \linkS4class{UnaryOperator} object.
#' @param values A list of arguments to the atom.
#' @describeIn UnaryOperator Applies the unary operator to the value.
setMethod("to_numeric", "UnaryOperator", function(object, values) {
  object@op_func(values[[1]])
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

setMethod("initialize", "NegExpression", function(.Object, ...) {
  callNextMethod(.Object, ..., op_name = "-", op_func = function(x) { -x })
})

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

#'
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
BinaryOperator <- setClass("BinaryOperator", representation(lh_exp = "ConstValORExpr", rh_exp = "ConstValORExpr", op_name = "character"), 
                                             prototype(op_name = "BINARY_OP"), contains = "AffAtom")

setMethod("initialize", "BinaryOperator", function(.Object, ..., lh_exp, rh_exp, op_name = "BINARY_OP") {
  .Object@lh_exp = lh_exp
  .Object@rh_exp = rh_exp
  .Object@op_name = op_name
  callNextMethod(.Object, ..., atom_args = list(.Object@lh_exp, .Object@rh_exp))
})

#' @param x,object A \linkS4class{BinaryOperator} object.
setMethod("name", "BinaryOperator", function(x) {
  pretty_args <- list()
  for(a in x@args) {
    if(is(a, "AddExpression") || is(a, "DivExpression"))
      pretty_args <- list(pretty_args, paste("(", name(a), ")", sep = ""))
    else
      pretty_args <- list(pretty_args, name(a))
  }
  paste(pretty_args[[1]], x@op_name, pretty_args[[2]])
})

#' @param values A list of arguments to the atom.
#' @describeIn BinaryOperator Apply the binary operator to the values.
setMethod("to_numeric", "BinaryOperator", function(object, values) {
  # values <- lapply(values, intf_convert_if_scalar)
  Reduce(object@op_name, values)
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

setMethod("initialize", "MulExpression", function(.Object, ...) {
  callNextMethod(.Object, ..., op_name = "%*%")
})

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

setMethod(".grad", "MulExpression", function(object, values) {
  if(is_constant(object@args[[1]]) || is_constant(object@args[[2]]))
    return(callNextMethod(object, values))
  
  # TODO: Verify that the following code is correct for non-affine arguments.
  X <- values[[1]]
  Y <- values[[2]]
  
  DX_rows <- size(object@args[[1]])
  cols <- size(object@args[[1]])
  
  # DX = [diag(Y11), diag(Y12), ...]
  #      [diag(Y21), diag(Y22), ...]
  #      [   ...       ...      ...]
  DX <- sparseMatrix(i = c(), j = c(), dims = c(DX_rows, cols))
  step <- dim(object@args[[1]])[1]
  for(k in 1:step)
    DX[seq(k, DX_rows, step), seq(k, cols, step)] <- Y
  if(length(dim(object@args[[2]])) == 1)
    cols <- 1
  else
    cols <- ncol(object@args[[2]])
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
#' The DivExpression class.
#'
#' This class represents one expression divided by another expression.
#'
#' @name DivExpression-class
#' @aliases DivExpression
#' @rdname DivExpression-class
.DivExpression <- setClass("DivExpression", contains = "BinaryOperator")
DivExpression <- function(lh_exp, rh_exp) { .DivExpression(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("initialize", "DivExpression", function(.Object, ...) {
  callNextMethod(.Object, ..., op_name = "/")
})

#' @param object A \linkS4class{DivExpression} object.
#' @describeIn DivExpression Is the left-hand expression quadratic and the right-hand expression constant?
setMethod("is_quadratic", "DivExpression", function(object) {
  is_quadratic(object@args[[1]]) && is_constant(object@args[[2]])
})

#' @describeIn DivExpression Is the expression quadratic of piecewise affine?
setMethod("is_qpwa", "DivExpression", function(object) {
  is_qpwa(object@args[[1]]) && is_constant(object@args[[2]])
})

#' @describeIn DivExpression The (row, col) dimensions of the left-hand expression.
setMethod("dim_from_args", "DivExpression", function(object) { dim(object@args[[1]]) })

#' @describeIn DivExpression Division is convex (affine) in its arguments only if the denominator is constant.
setMethod("is_atom_convex", "DivExpression", function(object) { is_constant(object@args[[2]]) && is_scalar(object@args[[2]]) })

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

#' @describeIn Multiply Multiplies the values elementwise.
setMethod("to_numeric", "Multiply", function(object, values) { values[[1]] * values[[2]] })

#' @describeIn Multiply The sum of the argument dimensions - 1.
setMethod("dim_from_args", "Multiply", function(object) { sum_dims(lapply(object@args, dim)) })

#' @describeIn Multiply Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Multiply", function(object) { TRUE })

#' @describeIn Multiply Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Multiply", function(object) { TRUE })

#' @describeIn Multiply Is the expression a positive semidefinite matrix?
setMethod("is_psd", "Multiply", function(object) {
  (is_psd(object@args[[1]]) && is_nsd(object@args[[2]])) || (is_nsd(object@args[[1]]) && is_psd(object@args[[2]]))
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
  apply(values[[1]], object@axis, base::cumsum)
})

#' @describeIn CumSum The dimensions of the atom.
setMethod("dim_from_args", "CumSum", function(object) { dim(object@args[[1]]) })

setMethod(".grad", "CumSum", function(object, values) {
  # TODO: This is inefficient
  val_dim <- dim(object@values[[1]])
  collapse <- setdiff(1:length(val_dim), object@axis)
  dim <- val_dim[collapse]
  mat <- matrix(0, nrow = dim, ncol = dim)
  mat[lower.tri(mat, diag = TRUE)] <- 1

  var <- Variable(dim(object@args[[1]]))
  if(object@axis == 2)
    grad <- .grad(new("MulExpression", lh_exp = mat, rh_exp = var), values)[[2]]
  else
    grad <- .grad(new("RMulExpression", lh_exp = var, rh_exp = t(mat)), values)[[1]]
  list(grad)
})

#' @describeIn CumSum Returns the axis being summed.
setMethod("get_data", "CumSum", function(object) { list(object@axis) })

CumSum.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  # Implicit O(n) definition:
  # X = Y[:1,:] - Y[1:,:]
  Y <- create_var(dim)
  axis <- data[[1]]
  collapse <- setdiff(1:length(dim), axis)
  dim <- dim[collapse]
  diff_mat <- get_diff_mat(dim, axis)
  diff_mat <- create_const(diff_mat, c(dim, dim), sparse = TRUE)

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

Diag <- function(expr) {
  expr <- as.Constant(expr)
  if(is_vector(expr))
    return(DiagVec(Vec(expr)))
  else if(ndim(expr) == 2 && nrow(expr) == ncol(expr))
    return(DiagMat(expr = expr))
  else
    stop("Argument to Diag must be a vector or square matrix.")
}

Diff <- function(x, lag = 1, k = 1, axis = 1) {
  x <- as.Constant(x)
  if((axis == 1 && ndim(x) < 2) || ndim(x) == 0)
    stop("Invalid axis given input dimensions.")
  else if(axis == 2)
    x <- t(x)
  
  x_dim <- dim(x)
  collapse <- setdiff(1:length(x_dim), axis)
  m <- x_dim[collapse]
  if(k <= 0 || k >= m)
    stop("Must have k > 0 and x must have < k elements along collapsed axis.")
  if(lag <= 0 || lag >= m)
    stop("Must have lag > 0 and x must have < lag elements along collapsed axis.")

  d <- x
  len <- m
  for(i in 1:k) {
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
  for(idx in length(arg_list)) {
    arg <- arg_list[[idx]]
    if(ndim(arg) == 0)
      arg_list[[idx]] <- as.vector(arg)
  }
  .HStack(atom_args = arg_list)
}

#' @param object A \linkS4class{HStack} object.
#' @param values A list of arguments to the atom.
#' @describeIn HStack Horizontally concatenate the values using \code{cbind}.
setMethod("to_numeric", "HStack", function(object, values) { Reduce("cbind", values) })

#' @describeIn HStack The dimensions of the atom.
setMethod("dim_from_args", "HStack", function(object) {
  if(ndim(object@args[[1]]) == 1)
    # return(c(sum(sapply(object@args, size)), NA))
    return(c(sum(sapply(object@args, size)), 1))
  else{
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
  .Object@expr
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
.SpecialIndex <- setClass("SpecialIndex", representation(expr = "Expression", key = "list", .select_mat = "numeric", .dim = "numeric"),
                          prototype(.select_mat = NA_real_, .dim = NA_real_), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} representing a vector or matrix.
#' @param key A list containing the start index, end index, and step size of the slice.
#' @rdname SpecialIndex-class
SpecialIndex <- function(expr, key) { .SpecialIndex(expr = expr, key = key) }

setMethod("initialize", "SpecialIndex", function(.Object, ..., expr, key) {
  .Object@key <- key
  row <- key[[1]]
  col <- key[[2]]
  
  # Order the entries of expr and select them using key.
  expr_dim <- dim(expr)
  expr_size <- size(expr)
  
  idx_mat <- seq(expr_size)
  idx_mat <- matrix(idx_mat, nrow = expr_dim[1], ncol = expr_dim[2])
  if(is.matrix(row) && is.null(col))
    select_mat <- idx_mat[row]
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

setMethod("name", "SpecialIndex", function(x) { paste(name(x@args[[1]]), as.character(x@key)) })

#' @param object An \linkS4class{Index} object.
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
  if(is_scalar(Vec(v)) || is_scalar(as.Constant(idmat)))
    lowered <- Reshape(idmat * v, c(final_dim[1], final_dim[2]))
  else
    lowered <- Reshape(idmat %*% v, c(final_dim[1], final_dim[2]))
  return(grad(lowered))
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
  if(!is_constant(object@args[[1]]))
    stop("The first argument to Kron must be constant.")
  else if(ndim(object@args[[1]]) != 2 || ndim(object@args[[2]]) != 2)
    stop("Kron requires matrix arguments.")
})

#' @describeIn Kron The dimensions of the atom.
setMethod("dim_from_args", "Kron", function(object) {
  rows <- dim(object@args[[1]])[1] * dim(object@args[[2]])[1]
  cols <- dim(object@args[[1]])[2] * dim(object@args[[2]])[2]
  c(rows, cols)
})

#' @describeIn Kron The sign of the atom.
setMethod("sign_from_args", "Kron", function(object) { mul_sign(object@args[[1]], object@args[[2]]) })

#' @param idx An index into the atom.
#' @describeIn Kron Is the left-hand expression positive?
setMethod("is_incr", "Kron", function(object, idx) { is_nonneg(object@args[[1]]) })

#' @describeIn Kron Is the right-hand expression negative?
setMethod("is_decr", "Kron", function(object, idx) { is_nonpos(object@args[[2]]) })

Kron.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.kron(arg_objs[[1]], arg_objs[[2]], dim), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Kron The graph implementation of the atom.
setMethod("graph_implementation", "Kron", function(object, arg_objs, dim, data = NA_real_) {
  Kron.graph_implementation(arg_objs, dim, data)
})

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

#' @describeIn Promote Promotes the value to the new dimensions.
setMethod("to_numeric", "Promote", function(object, values) {
  array(1, dim = object@promoted_dim) * values[[1]]
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
Imag <- function(expr) { .Real(expr = expr) }

setMethod("initialize", "Real", function(.Object, ..., expr) {
  .Object@expr
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
.Reshape <- setClass("Reshape", representation(expr = "ConstValORExpr", new_dim = "numeric"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} or numeric matrix.
#' @param new_dim The new dimensions.
#' @rdname Reshape-class
Reshape <- function(expr, new_dim) { .Reshape(expr = expr, new_dim = new_dim) }

setMethod("initialize", "Reshape", function(.Object, ..., expr, new_dim) {
  .Object@new_dim <- new_dim
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{Reshape} object.
#' @param values A list of arguments to the atom.
#' @describeIn Reshape Reshape the value into the specified dimensions.
setMethod("to_numeric", "Reshape", function(object, values) {
  dim(values[[1]]) <- object@new_dim
  values[[1]]
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
setMethod("get_data", "Reshape", function(object) { list(object@new_dim) })

Reshape.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.reshape(arg_objs[[1]], dim), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Reshape The graph implementation of the atom.
setMethod("graph_implementation", "Reshape", function(object, arg_objs, dim, data = NA_real_) {
  Reshape.graph_implementation(arg_objs, dim, data)
})

#'
#' The SumEntries class.
#'
#' This class represents the sum of all entries in a vector or matrix.
#'
#' @slot expr An \linkS4class{Expression} representing a vector or matrix.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @name SumEntries-class
#' @aliases SumEntries
#' @rdname SumEntries-class
.SumEntries <- setClass("SumEntries", contains = c("AxisAtom", "AffAtom"))

#' @param expr An \linkS4class{Expression} representing a vector or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
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
  axis <- data[[1]]
  keepdims <- data[[2]]
  if(is.na(axis))
    obj <- lo.sum_entries(arg_objs[[1]], dim)
  else if(axis == 1) {
    if(keepdims)
      const_dim <- c(arg_objs[[1]]$dim[2], 1)
    else
      const_dim <- c(arg_objs[[1]]$dim[2], NA_integer_)
    ones <- create_const(array(1, dim = const_dim), const_dim)
    obj <- lo.rmul_expr(arg_objs[[1]], ones, dim)
  } else {   # axis == 2
    if(keepdims)
      const_dim <- c(1, arg_objs[[1]]$dim[1])
    else
      const_dim <- c(arg_objs[[1]]$dim[1], NA_integer_)
    ones <- create_const(array(1, dim = const_dim), const_dim)
    obj <- lo.mul_expr(ones, arg_objs[[1]], dim)
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
setMethod("dim_from_args", "Trace", function(object){ c() })

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
setMethod("to_numeric", "Transpose", function(object, values) { aperm(values[[1]], perm = object@axes) })

#' @describeIn Transpose Is the expression symmetric?
setMethod("is_symmetric", "Transpose", function(object) { is_symmetric(object@args[[1]]) })

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

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn UpperTri The graph implementation of the atom.
setMethod("graph_implementation", "UpperTri", function(object, arg_objs, dim, data = NA_real_) {
  UpperTri.graph_implementation(arg_objs, dim, data)
})

Vec <- function(X) {
  X <- as.Constant(X)
  Reshape(expr = X, new_dim = size(X))
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
setMethod("to_numeric", "VStack", function(object, values) { Reduce("rbind", values) })

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

setMethod("rbind2", signature(x = "Expression", y = "ANY"), function(x, y, ...) { VStack(x, y) })
setMethod("rbind2", signature(x = "ANY", y = "Expression"), function(x, y, ...) { VStack(x, y) })

Bmat <- function(block_lists) {
  row_blocks <- lapply(block_lists, function(blocks) { .HStack(atom_args = blocks) })
  .VStack(atom_args = row_blocks)
}
