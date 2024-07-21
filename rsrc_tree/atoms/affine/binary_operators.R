## CVXPY SOURCE: cvxpy/atoms/affine/binary_operators.py

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
