#'
#' The AffAtom class.
#'
#' This virtual class represents an affine atomic expression.
#'
#' @aliases AffAtom
#' @export
AffAtom <- setClass("AffAtom", contains = c("VIRTUAL", "Atom"))

setMethod("func_curvature", "AffAtom", function(object) { Curvature(curvature = CURV_AFFINE_KEY) })
setMethod("sign_from_args", "AffAtom", function(object) { 
  arg_signs <- lapply(object@.args, function(arg) { arg@dcp_attr@sign })
  Reduce("+", arg_signs)
})
setMethod("monotonicity", "AffAtom", function(object) { rep(INCREASING, length(object@.args)) })

#'
#' The AddExpression class.
#'
#' This class represents the sum of any number of expressions.
#'
#' @slot arg_groups A \code{list} of \S4class{Expression}s and numeric data.frame, matrix, or vector objects.
#' @aliases AddExpression
#' @export
AddExpression <- setClass("AddExpression", representation(arg_groups = "list"), prototype(arg_groups = list()), contains = "AffAtom")

setMethod("init_dcp_attr", "AddExpression", function(object) {
  arg_dcp <- lapply(object@.args, function(arg) { arg@dcp_attr })
  Reduce("+", arg_dcp)
})

setMethod("initialize", "AddExpression", function(.Object, ..., arg_groups = list()) {
  .Object@arg_groups <- arg_groups
  .Object <- callNextMethod(.Object, ..., .args = arg_groups)   # Casts R values to Constant objects
  .Object@.args <- lapply(.Object@.args, function(group) { if(is(group,"AddExpression")) group@.args else group })
  .Object@.args <- flatten_list(.Object@.args)   # Need to flatten list of expressions
  return(.Object)
})

AddExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  arg_objs <- lapply(arg_objs, function(arg) { if(size(arg) != size) promote(arg, size) else arg })
  list(sum_expr(arg_objs), list())
}

setMethod("graph_implementation", "AddExpression", function(object, arg_objs, size, data = NA_real_) {
  AddExpression.graph_implementation(arg_objs, size, data)
})

#'
#' The UnaryOperator class.
#'
#' This base class represents expressions involving unary operators.
#'
#' @slot expr The \S4class{Expression} that is being operated upon.
#' @slot op_name A \code{character} string indicating the unary operation.
#' @aliases UnaryOperator
#' @export
UnaryOperator <- setClass("UnaryOperator", representation(expr = "Expression", op_name = "character"), contains = "AffAtom")

setMethod("init_dcp_attr", "UnaryOperator", function(object) {
  .Primitive(object@op_name)(object@.args[[1]]@dcp_attr)
})

setMethod("initialize", "UnaryOperator", function(.Object, ..., expr, op_name) {
  .Object@expr = expr
  .Object@op_name = op_name
  callNextMethod(.Object, ..., .args = list(.Object@expr))
})

#'
#' The NegExpression class.
#'
#' This class represents the negation of an affine expression.
#'
#' @aliases NegExpression
#' @export
NegExpression <- setClass("NegExpression", contains = "UnaryOperator")

setMethod("initialize", "NegExpression", function(.Object, ...) {
  callNextMethod(.Object, ..., op_name = "-")
})

NegExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(neg_expr(arg_objs[[1]]), list())
}

setMethod("graph_implementation", "NegExpression", function(object, arg_objs, size, data = NA_real_) {
  NegExpression.graph_implementation(arg_objs, size, data)
})

#'
#' The BinaryOperator class.
#'
#' This base class represents expressions involving binary operators.
#'
#' @slot lh_exp The \S4class{Expression} on the left-hand side of the operator.
#' @slot rh_exp The \S4class{Expression} on the right-hand side of the operator.
#' @slot op_name A \code{character} string indicating the binary operation.
#' @aliases BinaryOperator
#' @export
BinaryOperator <- setClass("BinaryOperator", representation(lh_exp = "Expression", rh_exp = "Expression", op_name = "character"), contains = "AffAtom")

setMethod("init_dcp_attr", "BinaryOperator", function(object) {
  .Primitive(object@op_name)(object@.args[[1]]@dcp_attr, object@.args[[2]]@dcp_attr)
})

setMethod("validate_args", "BinaryOperator", function(object) {
  .Primitive(object@op_name)(object@.args[[1]]@dcp_attr@shape, object@.args[[2]]@dcp_attr@shape)
})

setMethod("initialize", "BinaryOperator", function(.Object, ..., lh_exp, rh_exp, op_name) {
  .Object@lh_exp = lh_exp
  .Object@rh_exp = rh_exp
  .Object@op_name = op_name
  callNextMethod(.Object, ..., .args = list(.Object@lh_exp, .Object@rh_exp))
})

#'
#' The MulExpression class.
#'
#' This class represents the product of two linear expressions.
#'
#' @aliases MulExpression
#' @export
MulExpression <- setClass("MulExpression", contains = "BinaryOperator")

setMethod("initialize", "MulExpression", function(.Object, ...) {
  callNextMethod(.Object, ..., op_name = "*")
})

MulExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  if(size[2] != 1 && all(size(arg_objs[[2]]) == c(1,1))) {
    arg <- promote(arg_objs[[2]], list(size[2], 1))
    arg_objs[2] <- diag_vec(arg)
  }
  list(mul_expr(arg_objs[[1]], arg_objs[[2]], size), list())
}

setMethod("graph_implementation", "MulExpression", function(object, arg_objs, size, data = NA_real_) {
  MulExpression.graph_implementation(arg_objs, size, data)
})

#'
#' The RMulExpression class.
#'
#' This class represents product of an expression with a constant on the right.
#'
#' @aliases RMulExpression
#' @export
RMulExpression <- setClass("RMulExpression", contains = "MulExpression")

RMulExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  if(size[1] != 1 && all(size(arg_objs[[1]]) == c(1,1))) {
    arg <- promote(arg_objs[[1]], list(size[1], 1))
    arg_objs[[1]] <- diag_vec(arg)
  }
  list(rmul_expr(arg_objs[[1]], arg_objs[[2]], size), list())
}

setMethod("graph_implementation", "RMulExpression", function(object, arg_objs, size, data = NA_real_) {
  RMulExpression.graph_implementation(arg_objs, size, data)
})

#'
#' The DivExpression class.
#'
#' This class represents one expression divided by another expression.
#'
#' @aliases DivExpression
#' @export
DivExpression <- setClass("DivExpression", contains = "BinaryOperator")

setMethod("initialize", "DivExpression", function(.Object, ...) {
  .Object@op_name = "/"
  callNextMethod(.Object, ..., op_name = "/")
})

DivExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(div_expr(arg_objs[[1]], arg_objs[[2]]), list())
}

setMethod("graph_implementation", "DivExpression", function(object, arg_objs, size, data = NA_real_) {
  DivExpression.graph_implementation(arg_objs, size, data)
})

#'
#' The Conv class.
#'
#' This class represents the 1-D discrete convolution of two vectors.
#'
#' @slot lh_exp An \S4class{Expression} or R numeric data representing the left-hand vector.
#' @slot rh_exp An \S4class{Expression} or R numeric data representing the right-hand vector.
#' @aliases Conv
#' @export
.Conv <- setClass("Conv", representation(lh_exp = "ConstValORExpr", rh_exp = "ConstValORExpr"), contains = "AffAtom")
Conv <- function(lh_exp, rh_exp) { .Conv(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("validate_args", "Conv", function(object) {
  if(!is_vector(object@.args[[1]]) || !is_vector(object@.args[[2]]))
    stop("The arguments to conv must resolve to vectors.")
  if(!is_constant(object@.args[[1]]))
    stop("The first argument to conv must be constant.")
})

setMethod("initialize", "Conv", function(.Object, ..., lh_exp, rh_exp) {
  .Object@lh_exp <- lh_exp
  .Object@rh_exp <- rh_exp
  callNextMethod(.Object, ..., .args = list(.Object@lh_exp, .Object@rh_exp))
})

setMethod("shape_from_args", "Conv", function(object) {
  lh_length <- size(object@.args[[1]])[1]
  rh_length <- size(object@.args[[2]])[1]
  Shape(rows = lh_length + rh_length - 1, cols = 1)
})

setMethod("sign_from_args", "Conv", function(object) {
  object@.args[[1]]@dcp_attr@sign * object@.args[[2]]@dcp_attr@sign
})

Conv.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(conv(arg_objs[[1]], arg_objs[[2]], size), list())
}

setMethod("graph_implementation", "Conv", function(object, arg_objs, size, data = NA_real_) {
  Conv.graph_implementation(arg_objs, size, data)
})

#'
#' The DiagVec class.
#'
#' This class represents the conversion of a vector into a diagonal matrix.
#'
#' @slot expr An \S4class{Expression} representing the vector to convert.
#' @aliases DiagVec
#' @export
.DiagVec <- setClass("DiagVec", representation(expr = "Expression"), contains = "AffAtom")
DiagVec <- function(expr) { .DiagVec(expr = expr) }

setMethod("initialize", "DiagVec", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., .args = list(.Object@expr))
})

setMethod("shape_from_args", "DiagVec", function(object) {
  rows <- size(object@.args[[1]])[1]
  Shape(rows = rows, cols = rows)
})

DiagVec.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(diag_vec(arg_objs[[1]]), list())
}

setMethod("graph_implementation", "DiagVec", function(object, arg_objs, size, data = NA_real_) {
  DiagVec.graph_implementation(arg_objs, size, data)
})

#'
#' The DiagMat class.
#'
#' This class represents the extraction of the diagonal from a square matrix.
#'
#' @slot expr An \S4class{Expression} representing the matrix whose diagonal we are interested in.
#' @aliases DiagMat
#' @export
.DiagMat <- setClass("DiagMat", representation(expr = "Expression"), contains = "AffAtom")
DiagMat <- function(expr) { .DiagMat(expr = expr) }

setMethod("initialize", "DiagMat", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., .args = list(.Object@expr))
})

setMethod("shape_from_args", "DiagMat", function(object) {
  rows <- size(object@.args[[1]])[1]
  Shape(rows = rows, cols = 1)
})

DiagMat.graph_implementation <- function(object, arg_objs, size, data = NA_real_) {
  list(diag_mat(arg_objs[[1]]), list())
}

setMethod("graph_implementation", "DiagMat", function(object, arg_objs, size, data = NA_real_) {
  DiagMat.graph_implementation(arg_objs, size, data)
})

#'
#' Matrix Diagonals
#' 
#' Extracts the diagonal from a matrix or makes a vector into a diagonal matrix.
#' 
#' @param x An \S4class{Expression} or numeric constant representing a vector or diagonal matrix.
setMethod("diag", "Expression", function(x) { 
  x <- as.Constant(x)
  if(is_vector(x)) {
    if(size(x)[2] == 1)
      return(DiagVec(x))
    else {
      x <- Reshape(x, size(x)[2], 1)
      return(DiagVec(x))
    }
  } else if(size(x)[1] == size(x)[2])
    return(DiagMat(x))
  else
    stop("Argument to diag must be a vector or square matrix")
})

# TODO: Overload diff function in R
Diff <- function(x, k = 1) {
  x <- as.Constant(x)
  m <- size(x)[1]
  n <- size(x)[2]
  
  if(k < 0 || k >= m || n != 1)
    stop("Must have k >= 0 and x must be a 1-D vector with < k elements")
  
  d <- x
  for(i in 1:k)
    d <- d[2:length(d)] - d[1:(length(d)-1)]
  d
}

.HStack <- setClass("HStack", contains = "AffAtom")

setMethod("validate_args", "HStack", function(object) {
  arg_cols <- sapply(object@.args, function(arg) { size(arg)[1] })
  if(max(arg_cols) != min(arg_cols))
    stop("All arguments to HStack must have the same number of rows")
})

setMethod("shape_from_args", "HStack", function(object) {
  arg_cols <- sapply(object@.args, function(arg) { size(arg)[1] })
  cols <- sum(arg_cols)
  rows <- size(object@.args[[1]])[1]
  Shape(rows = rows, cols = cols)
})

HStack.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(hstack(arg_objs, size), list())
}

setMethod("graph_implementation", "HStack", function(object, arg_objs, size, data = NA_real_) {
  HStack.graph_implementation(arg_objs, size, data)
})

.Index <- setClass("Index", representation(expr = "Expression", key = "numeric"), contains = "AffAtom")
Index <- function(expr, key) { .Index(expr = expr, key = key) }

setMethod("initialize", "Index", function(.Object, ..., expr, key) {
  .Object@key <- key   # TODO: Need to validate key
  .Object@expr <- expr
  callNextMethod(.Object, ..., .args = list(.Object@expr))
})

setMethod("shape_from_args", "Index", function(object) {
  size <- ku_size(object@key, object@.args[[1]]@dcp_attr@shape)
  Shape(size[1], size[2])
})

Index.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  obj <- index(arg_objs[[1]], size, data[[1]])
  list(obj, list())
}

setMethod("graph_implementation", "Index", function(object, arg_objs, size, data = NA_real_) {
  Index.graph_implementation(arg_objs, size, data)
})

Index.get_slice <- function(matrix, constraints, row_start, row_end, col_start, col_end) {
  key <- list(c(row_start, row_end, NA), c(col_start, col_end, NA))
  rows <- row_end - row_start
  cols <- col_end - col_start
  graph <- Index.graph_implementation(list(matrix), c(rows, cols), list(key))
  slc <- graph[[1]]
  idx_constr <- graph[[2]]
  constraints <- c(constraints, idx_constr)
  slc   # TODO: Should I return constraints as well?
}

Index.block_eq <- function(matrix, block, constraints, row_start, row_end, col_start, col_end) {
  key <- list(c(row_start, row_end, NA), c(col_start, col_end, NA))   # TODO: What is equivalent of Python slice object in R?
  rows <- row_end - row_start
  cols <- col_end - col_start
  if(!all(size(block) == c(rows, cols)))
    stop("Block must have rows = ", rows, " and cols = ", cols)
  graph <- Index.graph_implementation(list(matrix), c(rows, cols), list(key))
  slc <- graph[[1]]
  idx_constr <- graph[[2]]
  constraints <- c(constraints, create_leq(slc, block), idx_constr)
}

.Kron <- setClass("Kron", representation(lh_exp = "ConstValORExpr", rh_exp = "ConstValORExpr"), contains = "AffAtom")
Kron <- function(lh_exp, rh_exp) { .Kron(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("validate_args", "Kron", function(object) {
  if(!is_constant(object@.args[[1]]))
    stop("The first argument to Kron must be constant.")
})

setMethod("initialize", "Kron", function(.Object, ..., lh_exp, rh_exp) {
  .Object@lh_exp <- lh_exp
  .Object@rh_exp <- rh_exp
  callNextMethod(.Object, ..., .args = list(.Object@lh_exp, .Object@rh_exp))
})

setMethod("shape_from_args", "Kron", function(object) {
  rows <- size(object@.args[[1]])[1] * size(object@.args[[2]])[1]
  cols <- size(object@.args[[1]])[2] * size(object@.args[[2]])[2]
  Shape(rows = rows, cols = cols)
})

setMethod("sign_from_args", "Kron", function(object) {
  object@.args[[1]]@dcp_attr@sign * object@.args[[2]]@dcp_attr@sign
})

Kron.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(kron(arg_objs[[1]], arg_objs[[2]], size), list())
}

setMethod("graph_implementation", "Kron", function(object, arg_objs, size, data = NA_real_) {
  Kron.graph_implementation(arg_objs, size, data)
})

MulElemwise <- setClass("MulElemwise", representation(lh_const = "Expression", rh_exp = "Expression"), contains = "AffAtom")

setMethod("init_dcp_attr", "MulElemwise", function(object) {
  mul_elemwise(object@.args[[1]]@dcp_attr, object@.args[[2]]@dcp_attr)
})

setMethod("validate_args", "MulElemwise", function(object) {
  if(!is_constant(object@.args[[1]]))
    stop("The first argument to MulElemwise must be constant.")
})

setMethod("initialize", "MulElemwise", function(.Object, ..., lh_const, rh_exp) {
  .Object@lh_const <- lh_const
  .Object@rh_exp <- rh_exp
  callNextMethod(.Object, ..., .args = list(.Object@lh_const, .Object@rh_exp))
})

MulElemwise.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  if(size(arg_objs[[1]]) != size(arg_objs[[2]]))
    list(mul_expr(arg_objs[[1]], arg_objs[[2]], size), list())
  else
    list(mul_elemwise(arg_objs[[1]], arg_objs[[2]]), list())
}

setMethod("graph_implementation", "MulElemwise", function(object, arg_objs, size, data = NA_real_) {
  MulElemwise.graph_implementation(arg_objs, size, data)
})

#'
#' The Reshape class.
#'
#' This class represents the reshaping of an expression. The operator vectorizes the expression,
#' then unvectorizes it into the new shape. Entries are stored in column-major order.
#'
#' @slot expr The \S4class{Expression} to reshape.
#' @slot rows The new number of rows.
#' @slot cols The new number of columns.
#' @aliases Reshape
#' @export
.Reshape <- setClass("Reshape", representation(expr = "Expression", rows = "numeric", cols = "numeric"), contains = "AffAtom")
Reshape <- function(expr, rows, cols) { .Reshape(expr = expr, rows = rows, cols = cols) }

setMethod("validate_args", "Reshape", function(object) {
  old_len <- size(object@.args[[1]])[1] * size(object@.args[[1]])[2]
  new_len <- object@rows * object@cols
  if(old_len != new_len)
    stop(sprintf("Invalid reshape dimensions (%i, %i)", object@rows, object@cols))
})

setMethod("initialize", "Reshape", function(.Object, ..., expr, rows, cols) {
  .Object@rows <- rows
  .Object@cols <- cols
  .Object@expr <- expr
  callNextMethod(.Object, ..., .args = list(.Object@expr))
})

setMethod("shape_from_args", "Reshape", function(object) {
  Shape(rows = object@rows, cols = object@cols)
})

Reshape.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(reshape(arg_objs[[1]], size), list())
}

setMethod("graph_implementation", "Reshape", function(object, arg_objs, size, data = NA_real_) {
  Reshape.graph_implementation(arg_objs, size, data)
})

#'
#' The SumEntries class.
#'
#' This class represents sum of all the entries in an expression.
#'
#' @slot expr The \S4class{Expression} to sum the entries of.
#' @aliases SumEntries
#' @export
.SumEntries <- setClass("SumEntries", representation(expr = "ConstValORExpr"), contains = "AffAtom")
SumEntries <- function(expr) { .SumEntries(expr = expr) }
setMethod("initialize", "SumEntries", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., .args = list(.Object@expr))
})

setMethod("shape_from_args", "SumEntries", function(object){
  Shape(rows = 1, cols = 1)
})

SumEntries.grpah_implementation <- function(arg_objs, size, data = NA_real_) {
  list(sum_entries(arg_objs[[1]]), list())
}

setMethod("graph_implementation", "SumEntries", function(object, arg_objs, size, data = NA_real_) {
  SumEntries.graph_implementation(arg_objs, size, data)
})

#'
#' The Trace class.
#'
#' This class represents the sum of the diagonal entries in a matrix.
#'
#' @slot expr The \S4class{Expression} to sum the diagonal of.
#' @aliases Trace
#' @export
.Trace <- setClass("Trace", representation(expr = "Expression"), contains = "AffAtom")
Trace <- function(expr) { .Trace(expr = expr) }

setMethod("validate_args", "Trace", function(object) {
  size <- size(object@.args[[1]])
  if(size[1] != size[2])
    stop("Argument to trace must be a square matrix")
})

setMethod("initialize", "Trace", function(.Object, ..., expr) {
  .Object@expr = expr
  callNextMethod(.Object, ..., .args = list(.Object@expr))
})

setMethod("shape_from_args", "Trace", function(object){
  Shape(rows = 1, cols = 1)
})

Trace.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(trace(arg_objs[[1]]), list())
}

setMethod("graph_implementation", "Trace", function(object, arg_objs, size, data = NA_real_) {
  Trace.graph_implementation(arg_objs, size, data)
})

#'
#' The Transpose class.
#'
#' This class represents the matrix transpose.
#'
#' @aliases Transpose
#' @export
Transpose <- setClass("Transpose", contains = "AffAtom")

setMethod("shape_from_args", "Transpose", function(object) {
  obj_size = size(object@.args[[1]])
  Shape(rows = obj_size[2], cols = obj_size[1])
})

Transpose.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(transpose(arg_objs[[1]]), list())  
}

setMethod("graph_implementation", "Transpose", function(object, arg_objs, size, data = NA_real_) {
  Transpose.graph_implementation(arg_objs, size, data)
})

UpperTri <- setClass("UpperTri", representation(expr = "Expression"), contains = "AffAtom")

setMethod("validate_args", "UpperTri", function(object) {
  if(size(object@.args[[1]])[1] != size(object@.args[[2]])[2])
    stop("Argument to UpperTri must be a square matrix.")
})

setMethod("initialize", "UpperTri", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., .args = list(.Object@expr))
})

setMethod("shape_from_args", "UpperTri", function(object) {
  rows <- size(object@.args[[1]])[1]
  cols <- size(object@.args[[2]])[2]
  Shape(rows = floor(rows*(cols-1)/2), cols = 1)
})

UpperTri.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(upper_tri(arg_objs[[1]]), list())
}

setMethod("graph_implementation", "UpperTri", function(object, arg_objs, size, data = NA_real_) {
  UpperTri.graph_implementation(arg_objs, size, data)
})

Vec <- function(X) {
  X <- as.Constant(X)
  Reshape(expr = X, rows = size(X)[1] * size(X)[2], cols = 1)
}

VStack <- setClass("VStack", contains = "AffAtom")

setMethod("validate_args", "VStack", function(object) {
  arg_cols <- sapply(object@.args, function(arg) { size(arg)[2] })
  if(max(arg_cols) != min(arg_cols))
    stop("All arguments to VStack must have the same number of columns")
})

setMethod("shape_from_args", "VStack", function(object) {
  cols <- size(object@.args[[1]])[2]
  arg_rows <- sapply(object@.args, function(arg) { size(arg)[1] })
  rows <- sum(arg_rows)
  Shape(rows = rows, cols = cols)
})

VStack.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(vstack(arg_objs, size), list())  
}

setMethod("graph_implementation", "VStack", function(object, arg_objs, size, data = NA_real_) {
  VStack.graph_implementation(arg_objs, size, data)
})

