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
#' @describeIn AffAtom Affine atoms can be complex valued.
setMethod("allow_complex", "Atom", function(object) { TRUE })

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

#' @describeIn AffAtom Is every argument quadratic or piecewise affine?
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
      fake_args <- c(fake_args, list(create_var(shape(arg), idx)))
      var_offsets <- c(var_offsets, offset)
      var_names <- c(var_names, idx)
      offset <- offset + size(arg)
    }
  }
  names(var_offsets) <- var_names
  graph <- graph_implementation(object, fake_args, shape(object), get_data(object))
  fake_expr <- graph[[1]]

  # Get the matrix representation of the function.
  prob_mat <- get_problem_matrix(list(create_eq(fake_expr)), var_offsets, NA)
  V <- prob_mat[[1]]
  I <- prob_mat[[2]] + 1   # TODO: R uses 1-indexing, but get_problem_matrix returns with 0-indexing
  J <- prob_mat[[3]] + 1
  shape <- c(offset, size(object))
  stacked_grad <- sparseMatrix(i = J, j = I, x = V, dims = shape)

  # Break up into per argument matrices.
  grad_list <- list()
  start <- 1
  for(arg in object@args) {
    if(is_constant(arg)) {
      grad_shape <- c(size(arg), shape[2])
      if(all(grad_shape == c(1,1)))
        grad_list <- c(grad_list, list(0))
      else
        grad_list <- c(grad_list, list(sparseMatrix(i = c(), j = c(), dims = grad_shape)))
    } else {
      stop <- start + size(arg)
      if(stop == start)
        grad_list <- c(grad_list, list(sparseMatrix(i = c(), j = c(), dims = c(0, shape[2]))))
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
AddExpression <- setClass("AddExpression", representation(arg_groups = "list"), prototype(arg_groups = list()), contains = "AffAtom")

setMethod("initialize", "AddExpression", function(.Object, ..., arg_groups = list()) {
  .Object@arg_groups <- arg_groups
  .Object <- callNextMethod(.Object, ..., args = arg_groups)   # Casts R values to Constant objects
  .Object@args <- lapply(.Object@args, function(group) { if(is(group,"AddExpression")) group@args else group })
  .Object@args <- flatten_list(.Object@args)   # Need to flatten list of expressions
})

#' @param x,object An \linkS4class{AddExpression} object.
#' @describeIn AddExpression The shape of the expression.
setMethod("shape_from_args", "AddExpression", function(object) { sum_shapes(lapply(object@args, shape)) })

setMethod("name", "AddExpression", function(x) {
  paste(sapply(object@args, as.character), collapse = " + ")
})

#' @param values A list of arguments to the atom.
#' @describeIn AddExpression Sum all the values.
setMethod("to_numeric", "AddExpression", function(object, values) {
  # values <- lapply(values, intf_convert_if_scalar)
  Reduce("+", values)
})

setMethod("is_symmetric", "AddExpression", function(object) {
  symm_args <- all(sapply(object@args, is_symmetric))
  return(shape(object)[1] == shape(object)[2] && symm_args)
})

setMethod("is_hermitian", "AddExpression", function(object) {
  herm_args <- all(sapply(object@args, is_hermitian))
  return(shape(object)[1] == shape(object)[2] && herm_args)
})

AddExpression.graph_implementation <- function(arg_objs, shape, data = NA_real_) {
  arg_objs <- lapply(arg_objs, function(arg) {
    if(!all(arg$shape == shape) && lo.is_scalar(arg)) 
      lo.promote(arg, shape) 
    else 
      arg })
  list(lo.sum_expr(arg_objs), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param shape A vector with two elements representing the shape of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn AddExpression The graph implementation of the expression.
setMethod("graph_implementation", "AddExpression", function(object, arg_objs, shape, data = NA_real_) {
  AddExpression.graph_implementation(arg_objs, shape, data)
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
  callNextMethod(.Object, ..., args = list(expr))
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
NegExpression <- setClass("NegExpression", contains = "UnaryOperator")

setMethod("initialize", "NegExpression", function(.Object, ...) {
  callNextMethod(.Object, ..., op_name = "-", op_func = function(x) { -x })
})

#' @describeIn NegExpression The (row, col) shape of the expression.
setMethod("shape_from_args", "NegExpression", function(object) { shape(object@args[[1]]) })

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

NegExpression.graph_implementation <- function(arg_objs, shape, data = NA_real_) {
  list(lo.neg_expr(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param shape A vector with two elements representing the shape of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn NegExpression The graph implementation of the expression.
setMethod("graph_implementation", "NegExpression", function(object, arg_objs, shape, data = NA_real_) {
  NegExpression.graph_implementation(arg_objs, shape, data)
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
  callNextMethod(.Object, ..., args = list(.Object@lh_exp, .Object@rh_exp))
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
#' See \linkS4class{MulElemwise} for the elementwise product.
#'
#' @seealso \linkS4class{MulElemwise}
#' @name MulExpression-class
#' @aliases MulExpression
#' @rdname MulExpression-class
MulExpression <- setClass("MulExpression", contains = "BinaryOperator")

setMethod("initialize", "MulExpression", function(.Object, ...) {
  callNextMethod(.Object, ..., op_name = "%*%")
})

#' @param object A \linkS4class{MulExpression} object.
#' @param values A list of arguments to the atom.
#' @describeIn MulExpression Matrix multiplication.
setMethod("to_numeric", "MulExpression", function(object, values) {
  return(values[[1]] %*% values[[2]])
})

#' @describeIn MulExpression The (row, col) shape of the expression.
setMethod("shape_from_args", "MulExpression", function(object) { mul_shapes(shape(object@args[[1]]), shape(object@args[[2]])) })

#' @describeIn MulExpression Multiplication is convex (affine) in its arguments only if one of the arguments is constant.
setMethod("is_atom_convex", "MulExpression", function(object) {
  is_constant(object@args[[1]]) || is_constant(object@args[[2]])
})

#' @describeIn MulExpression If the multiplication atom is convex, then it is affine.
setMethod("is_atom_concave", "MulExpression", function(object) { is_atom_convex(object) })

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
  step <- shape(object@args[[1]])[1]
  for(k in 1:step) {
    DX[seq(k, DX_rows, step), seq(k, cols, step)] <- Y
  if(length(shape(object@args[[2]])) == 1)
    cols <- 1
  else
    cols <- shape(object@args[[2]])[2]
  DY <- Matrix(bdiag(lapply(1:cols, function(k) { t(X) })), sparse = TRUE)
  return(list(DX, DY))
})

MulExpression.graph_implementation <- function(arg_objs, shape, data = NA_real_) {
  # Promote the right-hand side to a diagonal matrix if necessary
  lhs <- arg_objs[[1]]
  rhs <- arg_objs[[2]]
  if(lo.is_const(lhs))
    return(list(lo.mul_expr(lhs, rhs, shape), list()))
  else if(lo.is_const(rhs))
    return(list(lo.rmul_expr(lhs, rhs, shape), list()))
  else
    stop("Product of two non-constant expressions is not DCP.")
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param shape A vector with two elements representing the shape of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn MulExpression The graph implementation of the expression.
setMethod("graph_implementation", "MulExpression", function(object, arg_objs, shape, data = NA_real_) {
  MulExpression.graph_implementation(arg_objs, shape, data)
})

#'
#' The DivExpression class.
#'
#' This class represents one expression divided by another expression.
#'
#' @name DivExpression-class
#' @aliases DivExpression
#' @rdname DivExpression-class
DivExpression <- setClass("DivExpression", contains = "BinaryOperator")

setMethod("initialize", "DivExpression", function(.Object, ...) {
  callNextMethod(.Object, ..., op_name = "/")
})

#' @param object A \linkS4class{DivExpression} object.
#' @describeIn DivExpression Is the left-hand expression quadratic and the right-hand expression constant?
setMethod("is_quadratic", "DivExpression", function(object) {
  is_quadratic(object@args[[1]]) && is_constant(object@args[[2]])
})

#' @describeIn DivExpression Is the expression quadratic or piecewise affine?
setMethod("is_qpwa", "DivExpression", function(object) {
  is_qpwa(object@args[[1]]) && is_constant(object@args[[2]])
})

#' @describeIn DivExpression The (row, col) shape of the left-hand expression.
setMethod("shape_from_args", "DivExpression", function(object) { shape(object@args[[1]]) })

#' @param idx An index into the atom.
#' @describeIn DivExpression Is the right-hand expression positive?
setMethod("is_incr", "DivExpression", function(object, idx) { is_nonneg(object@args[[2]]) })

#' @describeIn DivExpression Is the right-hand expression negative?
setMethod("is_decr", "DivExpression", function(object, idx) { is_nonpos(object@args[[2]]) })

DivExpression.graph_implementation <- function(arg_objs, shape, data = NA_real_) {
  list(lo.div_expr(arg_objs[[1]], arg_objs[[2]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param shape A vector with two elements representing the shape of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn DivExpression The graph implementation of the expression.
setMethod("graph_implementation", "DivExpression", function(object, arg_objs, shape, data = NA_real_) {
  DivExpression.graph_implementation(arg_objs, shape, data)
})

# Multiplies two expressions elementwise.
.Multiply <- setClass("Multiply", contains = "MulExpression")
Multiply <- function(lh_exp, rh_exp) { .Multiply(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("initialize", "Multiply", function(.Object, ..., lh_exp, rh_exp) {
  lh_exp <- multiply.cast_to_const(lh_exp)
  rh_exp <- multiply.cast_to_const(rh_exp)
  if(is_scalar(lh_exp) && !is_scalar(rh_exp))
    lh_exp <- promote(lh_exp, shape(rh_exp))
  else if(is_scalar(rh_exp) && !is_scalar(lh_exp))
    rh_exp <- promote(rh_exp, shape(lh_exp))
  callNextMethod(.Object, ..., lh_exp = lh_exp, rh_exp = rh_exp)
})

#' @describeIn Multiply Multiples the values elementwise.
setMethod("to_numeric", "Multiply", function(object, values) { values[[1]] * values[[2]] })

#' @describeIn Multiply The sum of the argument dimensions - 1.
setMethod("shape_from_args", "Multiply", function(object) { sum_shapes(lapply(object@args, shape)) })

Multiply.graph_implementation <- function(arg_objs, shape, data = NA_real_) {
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
#' @param shape A vector with two elements representing the shape of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Multiply The graph implementation of the expression.
setMethod("graph_implementation", "Multiply", function(object, arg_objs, shape, data = NA_real_) {
  Multiply.graph_implementation(arg_objs, shape, data)
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
  callNextMethod(.Object, ..., args = list(.Object@expr))
})

#' @param object A \linkS4class{Conjugate} object.
#' @param values A list of arguments to the atom.
#' @describeIn Conjugate Elementwise complex conjugate of the constant.
setMethod("to_numeric", "Conjugate", function(object, values) { Conj(values[[1]]) })

#' @describeIn Conjugate The (row, col) shape of the expression.
setMethod("shape_from_args", "Conjugate", function(object) { shape(object@args[[1]]) })

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
  callNextMethod(.Object, ..., args = list(.Object@lh_exp, .Object@rh_exp))
})

#' @describeIn Conv Check both arguments are vectors and the first is a constant.
setMethod("validate_args", "Conv", function(object) {
  if(!is_vector(object@args[[1]]) || !is_vector(object@args[[2]]))
    stop("The arguments to Conv must resolve to vectors.")
  if(!is_constant(object@args[[1]]))
    stop("The first argument to Conv must be constant.")
})

#' @param object A \linkS4class{Conv} object.
#' @param values A list of arguments to the atom.
#' @describeIn Conv The convolution of the two values.
setMethod("to_numeric", "Conv", function(object, values) {
    .Call('_CVXR_cpp_convolve', PACKAGE = 'CVXR', as.vector(values[[1]]), as.vector(values[[2]]))
})

#' @describeIn Conv The size of the atom.
setMethod("size_from_args", "Conv", function(object) {
  lh_length <- size(object@args[[1]])[1]
  rh_length <- size(object@args[[2]])[1]
  c(lh_length + rh_length - 1, 1)
})

#' @describeIn Conv The sign of the atom.
setMethod("sign_from_args", "Conv", function(object) { mul_sign(object@args[[1]], object@args[[2]]) })

#' @param idx An index into the atom.
#' @describeIn Conv Is the left-hand expression positive?
setMethod("is_incr", "Conv", function(object, idx) { is_positive(object@args[[1]]) })

#' @describeIn Conv Is the left-hand expression negative?
setMethod("is_decr", "Conv", function(object, idx) { is_negative(object@args[[1]]) })

Conv.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(lo.conv(arg_objs[[1]], arg_objs[[2]], size), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Conv The graph implementation of the atom.
setMethod("graph_implementation", "Conv", function(object, arg_objs, size, data = NA_real_) {
  Conv.graph_implementation(arg_objs, size, data)
})

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
.CumSum <- setClass("CumSum", representation(axis = "numeric"), prototype(axis = 2),
                    validity = function(object) {
                      if(!(length(object@axis) == 1 && object@axis %in% c(1,2)))
                        stop("[CumSum: axis] axis must equal 1 (row) or 2 (column)")
                      return(TRUE)
                    }, contains = c("AxisAtom", "AffAtom"))

#' @param expr An \linkS4class{Expression} to be summed.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, and \code{2} indicates columns. The default is \code{2}.
#' @rdname CumSum-class
CumSum <- function(expr, axis = 2) { .CumSum(expr = expr, axis = axis) }

#' @param object A \linkS4class{CumSum} object.
#' @param values A list of arguments to the atom.
#' @describeIn CumSum The cumulative sum of the values along the specified axis.
setMethod("to_numeric", "CumSum", function(object, values) {
  if(object@axis == 1)
    t(apply(values[[1]], 1, base::cumsum))
  else    # axis = 2
    apply(values[[1]], 2, base::cumsum)
})

#' @describeIn CumSum The size of the atom.
setMethod("size_from_args", "CumSum", function(object) { size(object@args[[1]]) })

#' @describeIn CumSum Check that axis is either 1 or 2.
setMethod("validate_args", "CumSum", function(object) {
  if(length(object@axis) != 1 || !(object@axis %in% c(1,2)))
    stop("Invalid argument for axis: must equal 1 (row) or 2 (column)")
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

setMethod(".grad", "CumSum", function(object, values) {
  # TODO: This is inefficient
  if(object@axis == 1)
    dim <- dim(values[[1]])[2]
  else if(object@axis == 2)
    dim <- dim(values[[1]])[1]
  else
    stop("Invalid axis ", object@axis)
  mat <- matrix(0, nrow = dim, ncol = dim)
  mat[lower.tri(mat, diag = TRUE)] <- 1

  size <- size(object@args[[1]])
  var <- Variable(size[1], size[2])
  if(object@axis == 2)
    grad <- .grad(new("MulExpression", lh_exp = mat, rh_exp = var), values)[[2]]
  else
    grad <- .grad(new("RMulExpression", lh_exp = var, rh_exp = t(mat)), values)[[1]]
  list(grad)
})

CumSum.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  # Implicit 0(n) definition:
  # X = Y[:1,:] - Y[1:,:]
  Y <- create_var(size)
  axis <- data[[1]]
  if(axis == 2)
    dim <- size[1]
  else
    dim <- size[2]
  diff_mat <- get_diff_mat(dim, axis)
  diff_mat <- create_const(diff_mat, c(dim, dim), sparse = TRUE)

  if(axis == 2)
    diff <- lo.mul_expr(diff_mat, Y, size)
  else
    diff <- lo.rmul_expr(Y, diff_mat, size)
  list(Y, list(create_eq(arg_objs[[1]], diff)))
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn CumSum The graph implementation of the atom.
setMethod("graph_implementation", "CumSum", function(object, arg_objs, size, data = NA_real_) {
  CumSum.graph_implementation(arg_objs, size, data)
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
  callNextMethod(.Object, ..., args = list(.Object@expr))
})

#' @param object A \linkS4class{DiagVec} object.
#' @param values A list of arguments to the atom.
#' @describeIn DiagVec Convert the vector constant into a diagonal matrix.
setMethod("to_numeric", "DiagVec", function(object, values) { diag(as.vector(values[[1]])) })

#' @describeIn DiagVec The size of the atom.
setMethod("size_from_args", "DiagVec", function(object) {
  rows <- size(object@args[[1]])[1]
  c(rows, rows)
})

DiagVec.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(lo.diag_vec(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn DiagVec The graph implementation of the atom.
setMethod("graph_implementation", "DiagVec", function(object, arg_objs, size, data = NA_real_) {
  DiagVec.graph_implementation(arg_objs, size, data)
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
  callNextMethod(.Object, ..., args = list(.Object@expr))
})

#' @param object A \linkS4class{DiagMat} object.
#' @param values A list of arguments to the atom.
#' @describeIn DiagMat Extract the diagonal from a square matrix constant.
setMethod("to_numeric", "DiagMat", function(object, values) { diag(values[[1]]) })

#' @describeIn DiagMat The size of the atom.
setMethod("size_from_args", "DiagMat", function(object) {
  rows <- size(object@args[[1]])[1]
  c(rows, 1)
})

DiagMat.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(lo.diag_mat(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn DiagMat The graph implementation of the atom.
setMethod("graph_implementation", "DiagMat", function(object, arg_objs, size, data = NA_real_) {
  DiagMat.graph_implementation(arg_objs, size, data)
})

Diag <- function(expr) {
  expr <- as.Constant(expr)
  if(is_vector(expr)) {
    if(size(expr)[2] == 1)
      return(DiagVec(expr = expr))
    else {   # Convert a row vector to a column vector
      expr <- Reshape(expr, size(expr)[2], 1)
      return(DiagVec(expr = expr))
    }
  } else if(size(expr)[1] == size(expr)[2])
    return(DiagMat(expr = expr))
  else
    stop("Argument to Diag must be a vector or square matrix")
}

Diff <- function(x, lag = 1, k = 1, axis = 1) {
  x <- as.Constant(x)
  if(axis == 2)
    x <- t(x)

  size <- size(x)
  m <- size[1]
  n <- size[2]
  if(n != 1)
    stop("x must be a 1-D vector")
  if(k <= 0 || k >= m)
    stop("Must have 0 < k < number of elements in x")
  if(lag <= 0 || lag >= m)
    stop("Must have 0 < lag < number of elements in x")

  d <- x
  len <- m
  for(i in 1:k) {
    d <- d[(1+lag):m,] - d[1:(m-lag),]
    m <- m-1
  }

  if(axis == 2)
    t(d)
  else
    d
}

#'
#' The HStack class.
#'
#' Horizontal concatenation of values.
#'
#' @slot ... \linkS4class{Expression} objects or matrices. All arguments must have the same number of rows.
#' @name HStack-class
#' @aliases HStack
#' @rdname HStack-class
.HStack <- setClass("HStack", contains = "AffAtom")

#' @param ... \linkS4class{Expression} objects or matrices. All arguments must have the same number of rows.
#' @rdname HStack-class
HStack <- function(...) { .HStack(args = list(...)) }

#' @describeIn HStack Check all arguments have the same height.
setMethod("validate_args", "HStack", function(object) {
  arg_cols <- sapply(object@args, function(arg) { size(arg)[1] })
  if(max(arg_cols) != min(arg_cols))
    stop("All arguments to HStack must have the same number of rows")
})

#' @param object A \linkS4class{HStack} object.
#' @param values A list of arguments to the atom.
#' @describeIn HStack Horizontally concatenate the values using \code{cbind}.
setMethod("to_numeric", "HStack", function(object, values) { Reduce("cbind", values) })

#' @describeIn HStack The size of the atom.
setMethod("size_from_args", "HStack", function(object) {
  arg_cols <- sapply(object@args, function(arg) { size(arg)[2] })
  cols <- sum(arg_cols)
  rows <- size(object@args[[1]])[1]
  c(rows, cols)
})

HStack.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(lo.hstack(arg_objs, size), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn HStack The graph implementation of the atom.
setMethod("graph_implementation", "HStack", function(object, arg_objs, size, data = NA_real_) {
  HStack.graph_implementation(arg_objs, size, data)
})

setMethod("cbind2", signature(x = "Expression", y = "ANY"), function(x, y, ...) { HStack(x, y) })
setMethod("cbind2", signature(x = "ANY", y = "Expression"), function(x, y, ...) { HStack(x, y) })

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
  .Object@key <- ku_validate_key(key, size(expr))   # TODO: Double check key validation
  .Object@expr <- expr
  callNextMethod(.Object, ..., args = list(.Object@expr))
})

#' @param x,object An \linkS4class{Index} object.
#' @param values A list of arguments to the atom.
#' @describeIn Index The index/slice into the given value.
setMethod("to_numeric", "Index", function(object, values) {
  ku_slice_mat(values[[1]], object@key)
})

#' @describeIn Index The size of the atom.
setMethod("size_from_args", "Index", function(object) {
  ku_size(object@key, size(object@args[[1]]))
})

#' @describeIn Index A list containing \code{key}.
setMethod("get_data", "Index", function(object) { list(object@key) })

Index.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  obj <- lo.index(arg_objs[[1]], size, data[[1]])
  list(obj, list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Index The graph implementation of the atom.
setMethod("graph_implementation", "Index", function(object, arg_objs, size, data = NA_real_) {
  Index.graph_implementation(arg_objs, size, data)
})

#
# Get Special Slice Index
#
# Indexing using logical indexing or a list of indices.
#
# @param expr An \linkS4class{Expression} object.
# @param row The row index.
# @param col The column index.
# @return An \linkS4class{Expression} representing the index/slice.
# @rdname Index-get_special_slice
Index.get_special_slice <- function(expr, row, col) {
  expr <- as.Constant(expr)

  # Order the entries of expr and select them using key.
  expr_size <- size(expr)
  expr_prod <- prod(expr_size)

  idx_mat <- seq(expr_prod)
  idx_mat <- matrix(idx_mat, nrow = expr_size[1], ncol = expr_size[2])
  if(is.matrix(row) && is.null(col))
    select_mat <- idx_mat[row]
  else if(is.null(row) && !is.null(col))
    select_mat <- idx_mat[ , col, drop = FALSE]
  else if(!is.null(row) && is.null(col))
    select_mat <- idx_mat[row, , drop = FALSE]
  else
    select_mat <- idx_mat[row, col, drop = FALSE]

  if(!is.null(dim(select_mat)))
    final_size <- dim(select_mat)
  else   # Always cast 1-D arrays as column vectors
    final_size <- c(length(select_mat), 1)

  # Select the chosen entries from expr.
  select_vec <- as.vector(select_mat)
  ##select_vec <- as.matrix(select_mat, nrow=final_size[1L], ncol=final_size[2L])

  identity <- sparseMatrix(i = 1:expr_prod, j = 1:expr_prod, x = rep(1, expr_prod))
  idmat <- matrix(identity[select_vec, ], ncol = expr_prod)
  v <- Vec(expr)
  if(is_scalar(Vec(v)) || is_scalar(as.Constant(idmat)))
    Reshape(idmat * v, final_size[1], final_size[2])
  else
    Reshape(idmat %*% v, final_size[1], final_size[2])
}

#
# Get Index
#
# Returns a canonicalized index into a matrix.
#
# @param matrix A LinOp representing the matrix to be indexed.
# @param constraints A list of \linkS4class{Constraint} objects to append to.
# @param row The row index.
# @param col The column index.
# @return A list with the canonicalized index and updated constraints.
# @rdname Index-get_index
Index.get_index <- function(matrix, constraints, row, col) {
  key <- Key(row, col)
  graph <- Index.graph_implementation(list(matrix), c(1, 1), list(key))
  idx <- graph[[1]]
  idx_constr <- graph[[2]]
  constraints <- c(constraints, idx_constr)
  list(idx = idx, constraints = constraints)
}

#
# Index Block Equation
#
# Adds an equality setting a section of the matrix equal to a block. Assumes the block does not need to be promoted.
#
# @param matrix A LinOp representing the matrix to be indexed.
# @param block A LinOp representing the block in the block equality.
# @param constraints A list of \linkS4class{Constraint} objects to append to.
# @param row_start  The first row of the matrix section.
# @param row_end The last row of the matrix section.
# @param col_start The first column of the matrix section.
# @param col_end The last column of the matrix section.
# @return A list of \linkS4class{Constraint} objects.
# @rdname Index-block_eq
Index.block_eq <- function(matrix, block, constraints, row_start, row_end, col_start, col_end) {
  key <- Key(row_start:row_end, col_start:col_end)
  rows <- row_end - row_start + 1
  cols <- col_end - col_start + 1
  if(!all(size(block) == c(rows, cols)))
    stop("Block must have rows = ", rows, " and cols = ", cols)
  graph <- Index.graph_implementation(list(matrix), c(rows, cols), list(key))
  slc <- graph[[1]]
  idx_constr <- graph[[2]]
  constraints <- c(constraints, list(create_eq(slc, block)), idx_constr)
  constraints
}

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
  callNextMethod(.Object, ..., args = list(.Object@lh_exp, .Object@rh_exp))
})

#' @describeIn Kron Check both arguments are vectors and the first is a constant.
setMethod("validate_args", "Kron", function(object) {
  if(!is_constant(object@args[[1]]))
    stop("The first argument to Kron must be constant")
})

#' @param object A \linkS4class{Kron} object.
#' @param values A list of arguments to the atom.
#' @describeIn Kron The kronecker product of the two values.
setMethod("to_numeric", "Kron", function(object, values) {
  base::kronecker(values[[1]], values[[2]])
})

#' @describeIn Kron The size of the atom.
setMethod("size_from_args", "Kron", function(object) {
  rows <- size(object@args[[1]])[1] * size(object@args[[2]])[1]
  cols <- size(object@args[[1]])[2] * size(object@args[[2]])[2]
  c(rows, cols)
})

#' @describeIn Kron The sign of the atom.
setMethod("sign_from_args", "Kron", function(object) { mul_sign(object@args[[1]], object@args[[2]]) })

#' @param idx An index into the atom.
#' @describeIn Kron Is the left-hand expression positive?
setMethod("is_incr", "Kron", function(object, idx) { is_positive(object@args[[1]]) })

#' @describeIn Kron Is the right-hand expression negative?
setMethod("is_decr", "Kron", function(object, idx) { is_negative(object@args[[2]]) })

Kron.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(lo.kron(arg_objs[[1]], arg_objs[[2]], size), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Kron The graph implementation of the atom.
setMethod("graph_implementation", "Kron", function(object, arg_objs, size, data = NA_real_) {
  Kron.graph_implementation(arg_objs, size, data)
})

#'
#' The MulElemwise class.
#'
#' This class represents the elementwise multiplication of two expressions. The first expression must be constant.
#'
#' @slot lh_const A constant \linkS4class{Expression} or numeric value.
#' @slot rh_exp An \linkS4class{Expression}.
#' @name MulElemwise-class
#' @aliases MulElemwise
#' @rdname MulElemwise-class
.MulElemwise <- setClass("MulElemwise", representation(lh_const = "ConstValORExpr", rh_exp = "ConstValORExpr"), contains = "AffAtom")

#' @param lh_const A constant \linkS4class{Expression} or numeric value.
#' @param rh_exp An \linkS4class{Expression}.
#' @rdname MulElemwise-class
MulElemwise <- function(lh_const, rh_exp) { .MulElemwise(lh_const = lh_const, rh_exp = rh_exp) }

setMethod("initialize", "MulElemwise", function(.Object, ..., lh_const, rh_exp) {
  .Object@lh_const <- lh_const
  .Object@rh_exp <- rh_exp
  callNextMethod(.Object, ..., args = list(.Object@lh_const, .Object@rh_exp))
})

#' @describeIn MulElemwise Check the first argument is a constant.
setMethod("validate_args", "MulElemwise", function(object) {
  if(!is_constant(object@args[[1]]))
    stop("The first argument to MulElemwise must be constant.")
})

#' @param object A \linkS4class{MulElemwise} object.
#' @param values A list of arguments to the atom.
#' @describeIn MulElemwise Multiply the values elementwise.
setMethod("to_numeric", "MulElemwise", function(object, values) {
  values <- lapply(values, intf_convert_if_scalar)
  values[[1]] * values[[2]]
})

#' @describeIn MulElemwise The size of the atom.
setMethod("size_from_args", "MulElemwise", function(object) {
  sum_shapes(lapply(object@args, function(arg) { size(arg) }))
})

#' @describeIn MulElemwise The sign of the atom.
setMethod("sign_from_args", "MulElemwise", function(object) {
  mul_sign(object@args[[1]], object@args[[2]])
})

#' @param idx An index into the atom.
#' @describeIn MulElemwise Is the left-hand constant positive?
setMethod("is_incr", "MulElemwise", function(object, idx) { is_positive(object@args[[1]]) })

#' @describeIn MulElemwise Is the left-hand constant negative?
setMethod("is_decr", "MulElemwise", function(object, idx) { is_negative(object@args[[1]]) })

#' @describeIn MulElemwise Is the right-hand expression quadratic?
setMethod("is_quadratic", "MulElemwise", function(object) { is_quadratic(object@args[[2]]) })

MulElemwise.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  # One of the arguments is a scalar, so we can use normal multiplication
  if(any(arg_objs[[1]]$size != arg_objs[[2]]$size))
    MulExpression.graph_implementation(arg_objs, size, data)
  else
    list(lo.mul_elemwise(arg_objs[[1]], arg_objs[[2]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn MulElemwise The graph implementation of the atom.
setMethod("graph_implementation", "MulElemwise", function(object, arg_objs, size, data = NA_real_) {
  MulElemwise.graph_implementation(arg_objs, size, data)
})

#'
#' The Reshape class.
#'
#' This class represents the reshaping of an expression. The operator vectorizes the expression,
#' then unvectorizes it into the new shape. Entries are stored in column-major order.
#'
#' @slot expr An \linkS4class{Expression} or numeric matrix.
#' @slot rows The new number of rows.
#' @slot cols The new number of columns.
#' @name Reshape-class
#' @aliases Reshape
#' @rdname Reshape-class
.Reshape <- setClass("Reshape", representation(expr = "ConstValORExpr", rows = "numeric", cols = "numeric"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} or numeric matrix.
#' @param rows The new number of rows.
#' @param cols The new number of columns.
#' @rdname Reshape-class
Reshape <- function(expr, rows, cols) { .Reshape(expr = expr, rows = rows, cols = cols) }

setMethod("initialize", "Reshape", function(.Object, ..., expr, rows, cols) {
  .Object@rows <- rows
  .Object@cols <- cols
  .Object@expr <- expr
  callNextMethod(.Object, ..., args = list(.Object@expr))
})

#' @describeIn Reshape Check the new shape has the same number of entries as the old.
setMethod("validate_args", "Reshape", function(object) {
  old_len <- prod(size(object@args[[1]]))
  new_len <- object@rows * object@cols
  if(old_len != new_len)
    stop(sprintf("Invalid reshape dimensions (%i, %i)", object@rows, object@cols))
})

#' @param object A \linkS4class{Reshape} object.
#' @param values A list of arguments to the atom.
#' @describeIn Reshape Reshape the value into the specified dimensions.
setMethod("to_numeric", "Reshape", function(object, values) {
  dim(values[[1]]) <- c(object@rows, object@cols)
  values[[1]]
})

#' @describeIn Reshape The \code{c(rows, cols)} of the new expression.
setMethod("size_from_args", "Reshape", function(object) { c(object@rows, object@cols) })

#' @describeIn Reshape Returns \code{list(rows, cols)}.
setMethod("get_data", "Reshape", function(object) { list(object@rows, object@cols) })

Reshape.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(lo.reshape(arg_objs[[1]], size), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Reshape The graph implementation of the atom.
setMethod("graph_implementation", "Reshape", function(object, arg_objs, size, data = NA_real_) {
  Reshape.graph_implementation(arg_objs, size, data)
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
SumEntries <- function(expr, axis = NA_real_) { .SumEntries(expr = expr, axis = axis) }

#' @param object A \linkS4class{SumEntries} object.
#' @param values A list of arguments to the atom.
#' @describeIn SumEntries Sum the entries along the specified axis.
setMethod("to_numeric", "SumEntries", function(object, values) {
  if(is.na(object@axis))
    sum(values[[1]])
  else
    apply(values[[1]], object@axis, sum)
})

SumEntries.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  axis <- data[[1]]
  if(is.na(axis))
    obj <- lo.sum_entries(arg_objs[[1]])
  else if(axis == 1) {
    const_size <- c(arg_objs[[1]]$size[2], 1)
    ones <- create_const(matrix(1, nrow = const_size[1], ncol = const_size[2]), const_size)
    obj <- lo.rmul_expr(arg_objs[[1]], ones, size)
  } else {   # axis == 2
    const_size <- c(1, arg_objs[[1]]$size[1])
    ones <- create_const(matrix(1, nrow = const_size[1], ncol = const_size[2]), const_size)
    obj <- lo.mul_expr(ones, arg_objs[[1]], size)
  }
  list(obj, list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn SumEntries The graph implementation of the atom.
setMethod("graph_implementation", "SumEntries", function(object, arg_objs, size, data = NA_real_) {
  SumEntries.graph_implementation(arg_objs, size, data)
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
  callNextMethod(.Object, ..., args = list(.Object@expr))
})

#' @describeIn Trace Check the argument is a square matrix.
setMethod("validate_args", "Trace", function(object) {
  size <- size(object@args[[1]])
  if(size[1] != size[2])
    stop("Argument to trace must be a square matrix")
})

#' @param object A \linkS4class{Trace} object.
#' @param values A list of arguments to the atom.
#' @describeIn Trace Sum the diagonal entries.
setMethod("to_numeric", "Trace", function(object, values) { sum(diag(values[[1]])) })

#' @describeIn Trace The atom is a scalar.
setMethod("size_from_args", "Trace", function(object){ c(1, 1) })

Trace.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(lo.trace(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Trace The graph implementation of the atom.
setMethod("graph_implementation", "Trace", function(object, arg_objs, size, data = NA_real_) {
  Trace.graph_implementation(arg_objs, size, data)
})

#'
#' The Transpose class.
#'
#' This class represents the matrix transpose.
#'
#' @name Transpose-class
#' @aliases Transpose
#' @rdname Transpose-class
Transpose <- setClass("Transpose", contains = "AffAtom")

#' @param object A \linkS4class{Transpose} object.
#' @param values A list of arguments to the atom.
#' @describeIn Transpose The transpose of the given value.
setMethod("to_numeric", "Transpose", function(object, values) { t(values[[1]]) })

#' @describeIn Transpose The size of the atom.
setMethod("size_from_args", "Transpose", function(object) {
  size <- size(object@args[[1]])
  c(size[2], size[1])
})

Transpose.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(lo.transpose(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Transpose The graph implementation of the atom.
setMethod("graph_implementation", "Transpose", function(object, arg_objs, size, data = NA_real_) {
  Transpose.graph_implementation(arg_objs, size, data)
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
  callNextMethod(.Object, ..., args = list(.Object@expr))
})

#' @describeIn UpperTri Check the argument is a square matrix.
setMethod("validate_args", "UpperTri", function(object) {
  size <- size(object@args[[1]])
  if(size[1] != size[2])
    stop("Argument to UpperTri must be a square matrix.")
})

#' @param object An \linkS4class{UpperTri} object.
#' @param values A list of arguments to the atom.
#' @describeIn UpperTri Vectorize the upper triagonal entries.
setMethod("to_numeric", "UpperTri", function(object, values) {
  # Vectorize the upper triagonal entries
  tridx <- upper.tri(values[[1]], diag = FALSE)
  values[[1]][tridx]
})

#' @describeIn UpperTri The size of the atom.
setMethod("size_from_args", "UpperTri", function(object) {
  size <- size(object@args[[1]])
  rows <- size[1]
  cols <- size[2]
  c(rows*(cols-1)/2, 1)
})

UpperTri.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(lo.upper_tri(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn UpperTri The graph implementation of the atom.
setMethod("graph_implementation", "UpperTri", function(object, arg_objs, size, data = NA_real_) {
  UpperTri.graph_implementation(arg_objs, size, data)
})

Vec <- function(X) {
  X <- as.Constant(X)
  Reshape(expr = X, rows = prod(size(X)), cols = 1)
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
VStack <- function(...) { .VStack(args = list(...)) }

#' @describeIn VStack Check all arguments have the same width.
setMethod("validate_args", "VStack", function(object) {
  arg_cols <- sapply(object@args, function(arg) { size(arg)[2] })
  if(max(arg_cols) != min(arg_cols))
    stop("All arguments to VStack must have the same number of columns")
})

#' @param object A \linkS4class{VStack} object.
#' @param values A list of arguments to the atom.
#' @describeIn VStack Vertically concatenate the values using \code{rbind}.
setMethod("to_numeric", "VStack", function(object, values) { Reduce("rbind", values) })

#' @describeIn VStack The size of the atom.
setMethod("size_from_args", "VStack", function(object) {
  cols <- size(object@args[[1]])[2]
  arg_rows <- sapply(object@args, function(arg) { size(arg)[1] })
  rows <- sum(arg_rows)
  c(rows, cols)
})

VStack.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(lo.vstack(arg_objs, size), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn VStack The graph implementation of the atom.
setMethod("graph_implementation", "VStack", function(object, arg_objs, size, data = NA_real_) {
  VStack.graph_implementation(arg_objs, size, data)
})

setMethod("rbind2", signature(x = "Expression", y = "ANY"), function(x, y, ...) { VStack(x, y) })
setMethod("rbind2", signature(x = "ANY", y = "Expression"), function(x, y, ...) { VStack(x, y) })

Bmat <- function(block_lists) {
  row_blocks <- lapply(block_lists, function(blocks) { .HStack(args = blocks) })
  .VStack(args = row_blocks)
}
