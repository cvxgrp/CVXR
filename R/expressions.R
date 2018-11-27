#'
#' The Expression class.
#'
#' This class represents a mathematical expression.
#'
#' @name Expression-class
#' @aliases Expression
#' @rdname Expression-class
Expression <- setClass("Expression", contains = "Canonical")

setOldClass("data.frame")
setOldClass("matrix")
setOldClass("vector")
setClassUnion("ConstSparseVal", c("CsparseMatrix", "TsparseMatrix"))

setClassUnion("ConstVal", c("ConstSparseVal", "data.frame", "matrix", "vector", "numeric", "dMatrix"))
setClassUnion("ConstValORExpr", c("ConstVal", "Expression"))
setClassUnion("ListORExpr", c("list", "Expression"))
setClassUnion("ConstValListORExpr", c("ConstVal", "list", "Expression"))
setClassUnion("NumORgmp", c("numeric", "bigq", "bigz"))
setClassUnion("NumORNULL", c("numeric", "NULL"))

# Helper function since syntax is different for LinOp (list) vs. Expression object
#' @rdname size
setMethod("size", "ListORExpr", function(object) {
  if(is.list(object))
    object$size
  else
    size(object)
})

# Casts the second argument of a binary operator as an Expression
.cast_other <- function(binary_op) {
  cast_op <- function(object, other) {
    other <- as.Constant(other)
    binary_op(object, other)
  }
  cast_op
}

#' @param x,object An \linkS4class{Expression} object.
#' @describeIn Expression The value of the expression.
setMethod("value", "Expression", function(object) { stop("Unimplemented") })

#' @describeIn Expression The (sub/super)-gradient of the expression with respect to each variable.
setMethod("grad", "Expression", function(object) { stop("Unimplemented") })

#' @describeIn Expression A list of constraints describing the closure of the region where the expression is finite.
setMethod("domain", "Expression", function(object) { stop("Unimplemented") })

setMethod("show", "Expression", function(object) {
  cat("Expression(", curvature(object), ", ", sign(object), ", ", shape(object), ")", sep = "")
})

#' @rdname Expression-class
setMethod("as.character", "Expression", function(x) {
  paste("Expression(", curvature(x), ", ", sign(x), ", ", shape(x), ")", sep = "")
})

#' @describeIn Expression The string representation of the expression.
#' @export
setMethod("name", "Expression", function(object) { stop("Unimplemented") })

#' @describeIn Expression The expression itself.
setMethod("expr", "Expression", function(object) { object })

#'
#' Curvature of Expression
#'
#' The curvature of an expression.
#'
#' @param x An \linkS4class{Expression} object.
#' @return A string indicating the curvature of the expression, either "CONSTANT", "AFFINE", "CONVEX", "CONCAVE", or "UNKNOWN".
#' @docType methods
#' @rdname curvature
#' @export
setMethod("curvature", "Expression", function(object) {
  if(is_constant(object))
    curvature_str <- CONSTANT
  else if(is_affine(object))
    curvature_str <- AFFINE
  else if(is_convex(object))
    curvature_str <- CONVEX
  else if(is_concave(object))
    curvature_str <- CONCAVE
  else
    curvature_str <- UNKNOWN
  curvature_str
})

#' @describeIn Expression The expression is constant if it contains no variables or is identically zero.
setMethod("is_constant", "Expression", function(object) { length(variables(object)) == 0 || is_zero(object) || 0 %in% shape(object) })

#' @describeIn Expression The expression is affine if it is constant or both convex and concave.
setMethod("is_affine", "Expression", function(object) { is_constant(object) || (is_convex(object) && is_concave(object)) })

#' @describeIn Expression A logical value indicating whether the expression is convex.
setMethod("is_convex", "Expression", function(object) { stop("Unimplemented") })

#' @describeIn Expression A logical value indicating whether the expression is concave.
setMethod("is_concave", "Expression", function(object) { stop("Unimplemented") })

#' @describeIn Expression The expression is DCP if it is convex or concave.
setMethod("is_dcp", "Expression", function(object) { is_convex(object) || is_concave(object) })

#' @describeIn Expression A logical value indicating whether the expression is a Hermitian matrix.
setMethod("is_hermitian", "Expression", function(object) { is_real(object) && is_symmetric(object) })

#' @describeIn Expression A logical value indicating whether the expression is a positive semidefinite matrix.
setMethod("is_psd", "Expression", function(object) { FALSE })

#' @describeIn Expression A logical value indicating whether the expression is a negative semidefinite matrix.
setMethod("is_nsd", "Expression", function(object) { FALSE })

#' @describeIn Expression A logical value indicating whether the expression is quadratic.
setMethod("is_quadratic", "Expression", function(object) { FALSE })

#' @describeIn Expression A logical value indicating whether the expression is symmetric.
setMethod("is_symmetric", "Expression", function(object) { is_scalar(object) })

#' @describeIn Expression A logical value indicating whether the expression is piecewise linear.
setMethod("is_pwl", "Expression", function(object) { FALSE })

#' @describeIn Expression A logical value indicating whether the expression is quadratic or piecewise affine.
setMethod("is_qpwa", "Expression", function(object) { is_quadratic(object) || is_pwl(object) })

#'
#' Sign of Expression
#'
#' The sign of an expression.
#'
#' @param x An \linkS4class{Expression} object.
#' @return A string indicating the sign of the expression, either "ZERO", "NONNEG", "NONPOS", or "UNKNOWN".
#' @docType methods
#' @rdname sign
#' @export
setMethod("sign", "Expression", function(x) {
  if(is_zero(x))
    sign_str <- ZERO
  else if(is_nonneg(x))
    sign_str <- NONNEG
  else if(is_nonpos(x))
    sign_str <- NONPOS
  else
    sign_str <- UNKNOWN
  sign_str
})

#' @describeIn Expression The expression is zero if it is both nonnegative and nonpositive.
setMethod("is_zero", "Expression", function(object) { is_nonneg(object) && is_nonpos(object) })

#' @describeIn Expression A logical value indicating whether the expression is nonnegative.
setMethod("is_nonneg", "Expression", function(object) { stop("Unimplemented") })

#' @describeIn Expression A logical value indicating whether the expression is nonpositive.
setMethod("is_nonpos", "Expression", function(object) { stop("Unimplemented") })

#' @describeIn Expression The \code{c(row, col)} dimensions of the expression.
setMethod("shape", "Expression", function(object) { stop("Unimplemented") })

#' @describeIn Expression A logical value indicating whether the expression is real.
setMethod("is_real", "Expression", function(object) { !is_complex(object) })

#' @describeIn Expression A logical value indicating whether the expression is imaginary.
setMethod("is_imag", "Expression", function(object) { stop("Unimplemented") })

#' @describeIn Expression A logical value indicating whether the expression is complex.
setMethod("is_complex", "Expression", function(object) { stop("Unimplemented") })

#' @describeIn Expression The number of entries in the expression.
setMethod("size", "Expression", function(object) { as.integer(prod(shape(object))) })

#' @describeIn Expression The number of dimensions in the expression's shape.
setMethod("ndim", "Expression", function(object) { length(shape(object)) })

#' @describeIn Expression Vectorizes the expression.
setMethod("flatten", "Expression", function(object) { Vec(object) })

#' @describeIn Expression A logical value indicating whether the expression is a scalar.
setMethod("is_scalar", "Expression", function(object) { all(shape(object) == 1) })

#' @describeIn Expression A logical value indicating whether the expression is a row or column vector.
setMethod("is_vector", "Expression", function(object) { ndim(object) <= 1 || (ndim(object) == 1 && min(shape(object)) == 1) })

#' @describeIn Expression A logical value indicating whether the expression is a matrix.
setMethod("is_matrix", "Expression", function(object) { ndim(object) == 2 && shape(object)[1] > 1 && shape(object)[2] > 1 })

#' @describeIn Expression Number of rows in the expression.
#' @export
setMethod("nrow", "Expression", function(x) { shape(x)[1] })

#' @describeIn Expression Number of columns in the expression.
#' @export
setMethod("ncol", "Expression", function(x) { shape(x)[2] })

# Slice operators
#' @param i,j The row and column indices of the slice.
#' @param ... (Unimplemented) Optional arguments.
#' @param drop (Unimplemented) A logical value indicating whether the result should be coerced to the lowest possible dimension.
#' @rdname Index-class
#' @export
setMethod("[", signature(x = "Expression", i = "missing", j = "missing", drop = "ANY"), function(x, i, j, ..., drop) { x })

#' @rdname Index-class
#' @export
setMethod("[", signature(x = "Expression", i = "index", j = "missing", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  if(is_vector(x) && size(x)[1] < size(x)[2])
    Index.get_special_slice(x, NULL, i)   # If only first index given, apply it along longer dimension of vector
  else
    Index.get_special_slice(x, i, NULL)
})

#' @rdname Index-class
#' @export
setMethod("[", signature(x = "Expression", i = "missing", j = "index", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  Index.get_special_slice(x, NULL, j)
})

#' @rdname Index-class
#' @export
setMethod("[", signature(x = "Expression", i = "index", j = "index", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  Index.get_special_slice(x, i, j)
})

#' @rdname Index-class
#' @export
setMethod("[", signature(x = "Expression", i = "matrix", j = "index", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  Index.get_special_slice(x, i, j)
})

#' @rdname Index-class
#' @export
setMethod("[", signature(x = "Expression", i = "index", j = "matrix", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  Index.get_special_slice(x, i, j)
})

#' @rdname Index-class
#' @export
setMethod("[", signature(x = "Expression", i = "matrix", j = "matrix", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  Index.get_special_slice(x, i, j)
})

#' @rdname Index-class
#' @export
setMethod("[", signature(x = "Expression", i = "matrix", j = "missing", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  # This follows conventions in Matrix package, but differs from base handling of matrices
  Index.get_special_slice(x, i, NULL)
})

# #' @rdname Index-class
# setMethod("[", signature(x = "Expression", i = "ANY", j = "ANY", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
#  stop("Invalid or unimplemented Expression slice operation")
# })

# Arithmetic operators
#' @param e1,e2 The \linkS4class{Expression} objects or numeric constants to add.
#' @rdname AddExpression-class
setMethod("+", signature(e1 = "Expression", e2 = "missing"), function(e1, e2) { e1 })

#' @param e1,e2 The \linkS4class{Expression} objects or numeric constants to subtract.
#' @rdname NegExpression-class
setMethod("-", signature(e1 = "Expression", e2 = "missing"), function(e1, e2) { NegExpression(expr = e1) })

#' @rdname AddExpression-class
setMethod("+", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { AddExpression(arg_groups = list(e1, e2)) })

#' @rdname AddExpression-class
setMethod("+", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { AddExpression(arg_groups = list(e1, e2)) })

#' @rdname AddExpression-class
setMethod("+", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { e2 + e1 })

#' @rdname NegExpression-class
setMethod("-", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { e1 + NegExpression(expr = e2) })

#' @rdname NegExpression-class
setMethod("-", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { e1 + (-e2) })

#' @rdname NegExpression-class
setMethod("-", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { e1 + NegExpression(expr = e2) })

#' @param e1,e2 The \linkS4class{Expression} objects or numeric constants to multiply elementwise.
#' @docType methods
#' @rdname mul_elemwise
setMethod("*", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) {
  if(is_constant(e1))
    MulElemwise(lh_const = e1, rh_exp = e2)
  else if(is_constant(e2))
    MulElemwise(lh_const = e2, rh_exp = e1)
  else if(all(size(e1) == c(1,1)) && all(size(e2) == c(1,1)) && is_affine(e1) && is_affine(e2)) {
    warning("Forming a non-convex expression (affine) * (affine)")
    AffineProd(x = e1, y = e2)
  } else
    stop("Cannot multiply elementwise ", curvature(e1), " and ", curvature(e2))
})

#' @docType methods
#' @rdname mul_elemwise
setMethod("*", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { as.Constant(e2) * e1 })

#' @docType methods
#' @rdname mul_elemwise
setMethod("*", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { as.Constant(e1) * e2 })

#' @param e1,e2 The \linkS4class{Expression} objects or numeric constants to divide. The denominator, \code{e2}, must be a scalar constant.
#' @rdname DivExpression-class
setMethod("/", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) {
  if(is_constant(e2) && is_scalar(e2))
    DivExpression(lh_exp = e1, rh_exp = e2)
  else
    stop("Can only divide by a scalar constant")
})

#' @rdname DivExpression-class
setMethod("/", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { e1 / as.Constant(e2) })

#' @rdname DivExpression-class
setMethod("/", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { as.Constant(e1) / e2 })

#' @param e1 An \linkS4class{Expression} object to exponentiate.
#' @param e2 The power of the exponential. Must be a numeric scalar.
#' @docType methods
#' @rdname power
setMethod("^", signature(e1 = "Expression", e2 = "numeric"), function(e1, e2) {
  if(e2 == 2)
    Square(x = e1)
  else if(e2 == 0.5)
    Sqrt(x = e1)
  else
    Power(x = e1, p = e2)
})

# Matrix operators
#'
#' Matrix Transpose
#'
#' The transpose of a matrix.
#'
#' @param x An \linkS4class{Expression} representing a matrix.
#' @return An \linkS4class{Expression} representing the transposed matrix.
#' @docType methods
#' @aliases t
#' @rdname transpose
#' @method t Expression
#' @export
t.Expression <- function(x) { if(ndim(x) <= 1) x else Transpose(args = list(x)) }   # Need S3 method dispatch as well

#' @docType methods
#' @rdname transpose
#' @examples
#' x <- Variable(3, 4)
#' t(x)
#' @export
setMethod("t", signature(x = "Expression"), function(x) { if(ndim(x) <= 1) x else Transpose(args = list(x)) })

Conj.Expression <- function(z) { if(is_real(z)) z else Conjugate(z) }
setMethod("Conj", signature(z = "Expression"), function(z) { if(is_real(z)) z else Conjugate(z) })

#' @param x,y The \linkS4class{Expression} objects or numeric constants to multiply.
#' @rdname MulExpression-class
setMethod("%*%", signature(x = "Expression", y = "Expression"), function(x, y) {
  # if(is_scalar(x) || is_scalar(y))
  #  stop("Scalar operands are not allowed, use '*' instead")

  # Multiplying by a constant on the right is handled differently
  # from multiplying by a constant on the left
  if(is_constant(x)) {
    if(size(x)[1] == size(y)[1] && size(x)[2] != size(y)[1] && is(x, "Constant") && x@is_1D_array)
      x <- t(x)
    return(MulExpression(lh_exp = x, rh_exp = y))
  } else if(is_constant(y)) {
    # Having the constant on the left is more efficient
    if(is_scalar(x) || is_scalar(y))
      return(MulExpression(lh_exp = y, rh_exp = x))
    else
      return(RMulExpression(lh_exp = x, rh_exp = y))
  # When both expressions are not constant, allow affine * affine, but raise DCPError otherwise
  # Cannot multiply two non-constant expressions
  } else if(is_affine(x) && is_affine(y)) {
    warning("Forming a non-convex expression (affine) * (affine)")
    return(AffineProd(x = x, y = y))
  } else
    stop("Cannot multiply ", curvature(x), " and ", curvature(y))
})

#' @rdname MulExpression-class
setMethod("%*%", signature(x = "Expression", y = "ConstVal"), function(x, y) { x %*% as.Constant(y) })

#' @rdname MulExpression-class
setMethod("%*%", signature(x = "ConstVal", y = "Expression"), function(x, y) { as.Constant(x) %*% y })

# Comparison operators
#' @param e1,e2 The \linkS4class{Expression} objects or numeric constants to compare.
#' @rdname EqConstraint-class
setMethod("==", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { EqConstraint(lh_exp = e1, rh_exp = e2) })

#' @rdname EqConstraint-class
setMethod("==", signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 == as.Constant(e2) })

#' @rdname EqConstraint-class
setMethod("==", signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) == e2 })

#' @param e1,e2 The \linkS4class{Expression} objects or numeric constants to compare.
#' @rdname LeqConstraint-class
setMethod("<=", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { LeqConstraint(lh_exp = e1, rh_exp = e2) })

#' @rdname LeqConstraint-class
setMethod("<=", signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 <= as.Constant(e2) })

#' @rdname LeqConstraint-class
setMethod("<=", signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) <= e2 })

#' @rdname LeqConstraint-class
setMethod("<",  signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { stop("Unimplemented: Strict inequalities are not allowed.") })

#' @rdname LeqConstraint-class
setMethod("<",  signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 < as.Constant(e2) })

#' @rdname LeqConstraint-class
setMethod("<",  signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) < e2 })

#' @rdname LeqConstraint-class
setMethod(">=", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { e2 <= e1 })

#' @rdname LeqConstraint-class
setMethod(">=", signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 >= as.Constant(e2) })

#' @rdname LeqConstraint-class
setMethod(">=", signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) >= e2 })

#' @rdname LeqConstraint-class
setMethod(">",  signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { stop("Unimplemented: Strict inequalities are not allowed.") })

#' @rdname LeqConstraint-class
setMethod(">",  signature(e1 = "Expression", e2 = "ConstVal"),   function(e1, e2) { e1 > as.Constant(e2) })

#' @rdname LeqConstraint-class
setMethod(">",  signature(e1 = "ConstVal",   e2 = "Expression"), function(e1, e2) { as.Constant(e1) > e2 })

# Positive definite inequalities
#' @param e1,e2 The \linkS4class{Expression} objects or numeric constants to compare.
#' @docType methods
#' @rdname PSDConstraint-class
#' @export
setMethod("%>>%", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { PSDConstraint(e1, e2) })

#' @docType methods
#' @rdname PSDConstraint-class
#' @export
setMethod("%>>%", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { e1 %>>% as.Constant(e2) })

#' @docType methods
#' @rdname PSDConstraint-class
#' @export
setMethod("%>>%", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { as.Constant(e1) %>>% e2 })

#' @docType methods
#' @rdname PSDConstraint-class
#' @export
setMethod("%<<%", signature(e1 = "Expression", e2 = "Expression"), function(e1, e2) { PSDConstraint(e2, e1) })

#' @docType methods
#' @rdname PSDConstraint-class
#' @export
setMethod("%<<%", signature(e1 = "Expression", e2 = "ConstVal"), function(e1, e2) { e1 %<<% as.Constant(e2) })

#' @docType methods
#' @rdname PSDConstraint-class
#' @export
setMethod("%<<%", signature(e1 = "ConstVal", e2 = "Expression"), function(e1, e2) { as.Constant(e1) %<<% e2 })

#'
#' The Leaf class.
#'
#' This class represents a leaf node, i.e. a Variable, Constant, or Parameter.
#'
#' @slot args A list containing the arguments.
#' @name Leaf-class
#' @aliases Leaf
#' @rdname Leaf-class
Leaf <- setClass("Leaf", representation(shape = "numeric", value = "numeric", nonneg = "logical", nonpos = "logical",
                                        complex = "logical", imag = "logical", symmetric = "logical", diag = "logical",
                                        PSD = "logical", NSD = "logical", hermitian = "logical", boolean = "logical",
                                        integer = "logical", sparsity = "logical"),
                         prototype(value = NA_real_, nonneg = FALSE, nonpos = FALSE, complex = FALSE, imag = FALSE,
                                   symmetric = FALSE, diag = FALSE, PSD = FALSE, NSD = FALSE, hermitian = FALSE,
                                   boolean = FALSE, integer = FALSE, sparsity = NA), contains = "Expression")

setMethod("initialize", "Leaf", function(.Object, ..., shape, value = NA_real_, nonneg = FALSE, nonpos = FALSE, complex = FALSE, imag = FALSE, symmetric = FALSE, diag = FALSE, PSD = FALSE, NSD = FALSE, hermitian = FALSE, boolean = FALSE, integer = FALSE, sparsity = NA) {
  if(length(shape) > 2)
    stop("Expressions of dimension greater than 2 are not supported.")

  for(d in shape) {
    if(!is.integer(d) || d <= 0)
      stop("Invalid dimensions ", shape)
  }
  .Object@shape <- as.integer(shape)

  if((PSD || NSD || symmetric || diag || hermitian) && (length(shape) != 2 || shape[1] != shape[2]))
    stop("Invalid dimensions ", shape, ". Must be a square matrix.")

  # Process attributes.
  .Object@attributes <- list(nonneg = nonneg, nonpos = nonpos, complex = complex, imag = imag,
                             symmetric = symmetric, diag = diag, PSD = PSD, NSD = NSD, hermitian = hermitian,
                             boolean = as.logical(boolean), integer = integer, sparsity = sparsity)
  if(boolean) {
    if(!is.logical(boolean))
      .Object@boolean_idx <- boolean
    else
      .Object@boolean_idx <- do.call(expand.grid, lapply(shape, function(k) { 1:k }))
  } else
    .Object@boolean_idx <- list()

  if(integer) {
    if(!is.logical(integer))
      .Object@integer_idx <- integer
    else
      .Object@integer_idx <- do.call(expand.grid, lapply(shape, function(k) { 1:k }))
  } else
    .Object@integer_idx <- list()

  # Only one attribute can be TRUE (except boolean and integer).
  true_attr <- sum(unlist(.Object@attributes))
  if(boolean && integer)
    true_attr <- true_attr - 1
  if(true_attr > 1)
    stop("Cannot set more than one special attribute.")

  if(!is.na(value))
    .Object@value <- value
  .Object@args <- list()
  callNextMethod(.Object, ...)
})

setMethod("get_attr_str", "Leaf", function(object) {
  # Get a string representing the attributes
  attr_str <- ""
  for(attr in names(object@attributes)) {
    val <- object@attributes[attr]
    if(attr != "real" && !is.na(val))
      attr_str <- paste(attr_str, sprintf("%s=%s", attr, val), sep = ", ")
  }
  attr_str
})

#' @param object A \linkS4class{Leaf} object.
#' @describeIn Leaf Leaves are not copied.
setMethod("get_data", "Leaf", function(object) { })

#' @describeIn Leaf The dimensions of the leaf node.
setMethod("shape", "Leaf", function(object) { object@shape })

#' @describeIn Leaf List of \linkS4class{Variable} objects in the leaf node.
setMethod("variables", "Leaf", function(object) { list() })

#' @describeIn Leaf List of \linkS4class{Parameter} objects in the leaf node.
setMethod("parameters", "Leaf", function(object) { list() })

#' @describeIn Leaf List of \linkS4class{Constant} objects in the leaf node.
setMethod("constants", "Leaf", function(object) { list() })

#' @describeIn Leaf A logical value indicating whether the leaf node is convex.
setMethod("is_convex", "Leaf", function(object) { TRUE })

#' @describeIn Leaf A logical value indicating whether the leaf node is concave.
setMethod("is_concave", "Leaf", function(object) { TRUE })

#' @describeIn Leaf A logical value indicating whether the leaf node is nonnegative.
setMethod("is_nonneg", "Leaf", function(object) { object@attributes$nonneg || object@attributes$boolean })

#' @describeIn Leaf A logical value indicating whether the leaf node is nonpositive.
setMethod("is_nonpos", "Leaf", function(object) { object@attributes$nonpos })

#' @describeIn Leaf A logical value indicating whether the leaf node is hermitian.
setMethod("is_hermitian", "Leaf", function(object) {
  (is_real(object) && is_symmetric(object)) || object@attributes$hermitian || is_psd(object) || is_nsd(object)
})

#' @describeIn Leaf A logical value indicating whether the leaf node is symmetric.
setMethod("is_symmetric", "Leaf", function(object) {
  is_scalar(object) || any(sapply(c("diag", "symmetric", "PSD", "NSD"), function(key) { object@attributes[key] }))
})

#' @describeIn Leaf A logical value indicating whether the leaf node is imaginary.
setMethod("is_imag", "Leaf", function(object) { object@attributes$imag })

#' @describeIn Leaf A logical value indicating whether the leaf node is complex.
setMethod("is_complex", "Leaf", function(object) {
  object@attributes$complex || is_imag(object) || object@attributes$hermitian
})

#' @describeIn Leaf A list of constraints describing the closure of the region where the leaf node is finite. Default is the full domain.
setMethod("domain", "Leaf", function(object) {
  domain <- list()
  if(object@attributes$nonneg)
    domain <- c(domain, object >= 0)
  else if(object@attributes$nonpos)
    domain <- c(domain, object <= 0)
  else if(object@attributes$PSD)
    domain <- c(domain, object %>>% 0)
  else if(object@attributes$NSD)
    domain <- c(domain, object %<<% 0)
  return(domain)
})

#' @describeIn Leaf Project value onto the attribute set of the leaf.
setMethod("project", "Leaf", function(object, val) {
  if(!is_complex(object))
    val <- Re(val)

  if(object@attributes$nonpos && object@attributes$nonneg)
    return(0*val)
  else if(object@attributes$nonpos)
    return(pmin(val, 0))
  else if(object@attributes$nonneg)
    return(pmax(val, 0))
  else if(object@attributes$imag)
    return(Im(val))
  else if(object@attributes$complex)
    return(as.complex(val))
  else if(object@attributes$boolean)
    # TODO: Respect the boolean indices.
    return(round(pmax(pmin(val, 1), 0)))
  else if(object@attributes$integer)
    # TODO: Respect the integer indices. Also, variable may be integer in some indices and boolean in others.
    return(round(val))
  else if(object@attributes$diag) {
    val <- diag(val)
    return(sparseMatrix(i = 1:length(val), j = 1:length(val), x = val))
  } else if(object@attributes$hermitian)
    return(val + t(Conj(val))/2)
  else if(any(sapply(c("symmetric", "PSD", "NSD"), function(key) { object@attributes[key] }))) {
    val <- val + t(val)
    val <- val/2
    if(object@attributes$symmetric)
      return(val)

    wV <- eigen(val, symmetric = TRUE, only.values = FALSE)
    w <- wV$values
    V <- wV$vectors

    if(object@attributes$PSD) {
      bad <- w < 0
      if(!any(bad))
        return(val)
      w[bad] <- 0
    } else {   # NSD
      bad <- w > 0
      if(!any(bad))
        return(val)
      w[bad] <- 0
    }
    return((V %*% w) %*% t(V))
  } else
    return(val)
})

#' @param val The assigned value.
#' @describeIn Leaf Get the value of the leaf.
setMethod("value", "Leaf", function(object) { object@value })

#' @describeIn Leaf Set the value of the leaf.
setReplaceMethod("value", "Leaf", function(object, value) {
  object@value <- validate_val(object, value)
  return(object)
})

#' @describeIn Leaf Project and assign a value to the leaf.
setMethod("project_and_assign", "Leaf", function(object, val) {
  object@value <- project(object, val)
  return(object)
})

#' @describeIn Leaf Check that \code{val} satisfies symbolic attributes of leaf.
setMethod("validate_val", "Leaf", function(object, val) {
  if(!is.na(val)) {
    val <- intf_convert(val)
    if(any(intf_shape(val) != shape(object)))
      stop("Invalid dimensions ", intf_shape(val), " for value")
    projection <- project(object, val)
    delta <- abs(val - projection)

    if(is(delta, "sparseMatrix"))
      close_enough <- all(abs(delta@x) <= SPARSE_PROJECTION_TOL)
    else {
      delta <- as.matrix(delta)
      if(object@attributes$PSD || object@attributes$NSD)
        close_enough <- norm(delta, type = "2") <= PSD_NSD_PROJECTION_TOL
      else
        close_enough <- all(abs(delta) <= GENERAL_PROJECTION_TOL)
    }

    if(!close_enough) {
      if(object@attributes$nonneg)
        attr_str <- "nonnegative"
      else if(object@attributes$nonpos)
        attr_str <- "nonpositive"
      else if(object@attributes$diag)
        attr_str <- "diagonal"
      else if(object@attributes$PSD)
        attr_str <- "positive semidefinite"
      else if(object@attributes$NSD)
        attr_str <- "negative semidefinite"
      else if(object@attributes$imag)
        attr_str <- "imaginary"
      else {
        attr_str <- names(object@attributes)[unlist(object@attributes)]
        attr_str <- c(attr_str, "real")[1]
      }
      stop("Value must be ", attr_str)
    }
  }
  return(val)
})

#' @describeIn Leaf A logical value indicating whether the leaf node is a positive semidefinite matrix.
setMethod("is_psd", "Leaf", function(object) { object@attributes$PSD })

#' @describeIn Leaf A logical value indicating whether the leaf node is a negative semidefinite matrix.
setMethod("is_nsd", "Leaf", function(object) { object@attributes$NSD })

#' @describeIn Leaf Leaf nodes are always quadratic.
setMethod("is_quadratic", "Leaf", function(object) { TRUE })

#' @describeIn Leaf Leaf nodes are always piecewise linear.
setMethod("is_pwl", "Leaf", function(object) { TRUE })

setMethod("atoms", "Leaf", function(object) { list() })
