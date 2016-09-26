#'
#' The AffAtom class.
#'
#' This virtual class represents an affine atomic expression.
#'
#' @aliases AffAtom
#' @export
AffAtom <- setClass("AffAtom", contains = c("VIRTUAL", "Atom"))

setMethod("sign_from_args", "AffAtom", function(object) { sum_signs(object@args) })
setMethod("is_atom_convex", "AffAtom", function(object) { TRUE })
setMethod("is_atom_concave", "AffAtom", function(object) { TRUE })
setMethod("is_incr", "AffAtom", function(object, idx) { TRUE })
setMethod("is_decr", "AffAtom", function(object, idx) { FALSE })
setMethod("is_quadratic", "AffAtom", function(object) { all(sapply(object@args, function(arg) { is_quadratic(arg) })) })

.grad.AffAtom <- function(object, values) {
  # TODO: Should be a simple function in CVXcanon for this
  # Make a fake LinOp tree for the function
  fake_args <- list()
  var_offsets <- list()
  offset <- 0
  idx <- 1
  for(arg in object@args) {
    if(is_constant(arg))
      fake_args <- c(fake_args, list(create_const(value(arg), size(arg))))
    else {
      fake_args <- c(fake_args, list(create_var(size(arg), idx)))
      var_offsets[idx] <- offset
      offset <- offset + prod(size(arg))
    }
    idx <- idx + 1
  }
  graph <- graph_implementation(object, fake_args, size(object), get_data(object))
  fake_expr <- graph[[1]]
  
  # Get the matrix representation of the function
  prob_mat <- get_problem_matrix(list(create_eq(fake_expr)), var_offsets, NA)  # TODO: Call to canonInterface
  V <- prob_mat[[1]]
  I <- prob_mat[[2]]
  J <- prob_mat[[3]]
  shape <- c(offset, prod(size(object)))
  stacked_grad <- sparseMatrix(i = J, j = I, x = V, dims = shape)
  
  # Break up into per argument matrices
  grad_list <- list()
  start <- 1
  for(arg in object@args) {
    if(is_constant(arg)) {
      grad_shape <- c(prod(size(arg)), shape[2])
      if(all(grad_shape == c(1,1)))
        grad_list <- c(grad_list, list(0))
      else
        grad_list <- c(grad_list, list(Matrix(i = c(), j = c(), dims = grad_shape)))
    } else {
      stop <- start + prod(size(arg))
      grad_list <- c(grad_list, list(stacked_grad[start:stop,]))
      start <- stop
    }
  }
  grad_list
}

#'
#' The AddExpression class.
#'
#' This class represents the sum of any number of expressions.
#'
#' @slot arg_groups A \code{list} of \S4class{Expression}s and numeric data.frame, matrix, or vector objects.
#' @aliases AddExpression
#' @export
AddExpression <- setClass("AddExpression", representation(arg_groups = "list"), prototype(arg_groups = list()), contains = "AffAtom")

setMethod("initialize", "AddExpression", function(.Object, ..., arg_groups = list()) {
  .Object@arg_groups <- arg_groups
  .Object <- callNextMethod(.Object, ..., args = arg_groups)   # Casts R values to Constant objects
  .Object@args <- lapply(.Object@args, function(group) { if(is(group,"AddExpression")) group@args else group })
  .Object@args <- flatten_list(.Object@args)   # Need to flatten list of expressions
  return(.Object)
})

setMethod("to_numeric", "AddExpression", function(object, values) { Reduce("+", values)  })
setMethod("size_from_args", "AddExpression", function(object) { sum_shapes(lapply(object@args, function(arg) { size(arg) })) })

AddExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  arg_objs <- lapply(arg_objs, function(arg) { if(!all(arg$size == size)) promote(arg, size) else arg })
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

setMethod("initialize", "UnaryOperator", function(.Object, ..., expr, op_name) {
  .Object@expr <- expr
  .Object@op_name <- op_name
  callNextMethod(.Object, ..., args = list(.Object@expr))
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

setMethod("to_numeric", "NegExpression", function(object, values) { -values[[1]] })
setMethod("size_from_args", "NegExpression", function(object) { size(object@args[[1]]) })
setMethod("sign_from_args", "NegExpression", function(object) { c(is_negative(object@args[[1]]), is_positive(object@args[[1]])) })
setMethod("is_incr", "NegExpression", function(object, idx) { FALSE })
setMethod("is_decr", "NegExpression", function(object, idx) { TRUE })

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
BinaryOperator <- setClass("BinaryOperator", representation(lh_exp = "ConstValORExpr", rh_exp = "ConstValORExpr", op_name = "character"), contains = "AffAtom")

setMethod("initialize", "BinaryOperator", function(.Object, ..., lh_exp, rh_exp, op_name) {
  .Object@lh_exp = lh_exp
  .Object@rh_exp = rh_exp
  .Object@op_name = op_name
  callNextMethod(.Object, ..., args = list(.Object@lh_exp, .Object@rh_exp))
})

setMethod("to_numeric", "BinaryOperator", function(object, values) { Reduce(object@op_name, values) })
setMethod("sign_from_args", "BinaryOperator", function(object) { mul_sign(object@args[[1]], object@args[[2]]) })

#'
#' The MulExpression class.
#'
#' This class represents the matrix product of two linear expressions.
#' See MulElemwise for the elementwise product.
#'
#' @aliases MulExpression
#' @export
MulExpression <- setClass("MulExpression", contains = "BinaryOperator")

setMethod("initialize", "MulExpression", function(.Object, ...) {
  callNextMethod(.Object, ..., op_name = "%*%")
})

setMethod("validate_args", "MulExpression", function(object) {
  mul_shapes(size(object@args[[1]]), size(object@args[[2]]))
})

setMethod("size_from_args", "MulExpression", function(object) { mul_shapes(size(object@args[[1]]), size(object@args[[2]])) })
setMethod("is_incr", "MulExpression", function(object, idx) { is_positive(object@args[[1]]) })
setMethod("is_decr", "MulExpression", function(object, idx) { is_negative(object@args[[1]]) })

MulExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  # Promote the right-hand side to a diagonal matrix if necessary
  if(size[2] != 1 && all(arg_objs[[2]]$size == c(1,1))) {
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
#' This class represents the matrix product of an expression with a constant on the right.
#'
#' @aliases RMulExpression
#' @export
RMulExpression <- setClass("RMulExpression", contains = "MulExpression")

setMethod("is_incr", "RMulExpression", function(object, idx) { is_positive(object@args[[2]]) })
setMethod("is_decr", "RMulExpression", function(object, idx) { is_negative(object@args[[2]]) })

RMulExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  # Promote the left-hand side to a diagonal matrix if necessary
  if(size[1] != 1 && all(arg_objs[[1]]$size == c(1,1))) {
    arg <- promote(arg_objs[[1]], c(size[1], 1))
    arg_objs[1] <- diag_vec(arg)
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
  callNextMethod(.Object, ..., op_name = "/")
})

setMethod("is_quadratic", "DivExpression", function(object) {
  is_quadratic(object@args[[1]]) && is_constant(object@args[[2]])
})

setMethod("size_from_args", "DivExpression", function(object) { size(object@args[[1]]) })
setMethod("is_incr", "DivExpression", function(object, idx) { is_positive(object@args[[2]]) })
setMethod("is_decr", "DivExpression", function(object, idx) { is_negative(object@args[[2]]) })

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

setMethod("initialize", "Conv", function(.Object, ..., lh_exp, rh_exp) {
  .Object@lh_exp <- lh_exp
  .Object@rh_exp <- rh_exp
  callNextMethod(.Object, ..., args = list(.Object@lh_exp, .Object@rh_exp))
})

setMethod("validate_args", "Conv", function(object) {
  if(!is_vector(object@args[[1]]) || !is_vector(object@args[[2]]))
    stop("The arguments to Conv must resolve to vectors.")
  if(!is_constant(object@args[[1]]))
    stop("The first argument to Conv must be constant.")
})

setMethod("to_numeric", "Conv", function(object, values) {
  convolve(as.vector(values[[1]]), as.vector(values[[2]]))
})

setMethod("size_from_args", "Conv", function(object) {
  lh_length <- size(object@args[[1]])[1]
  rh_length <- size(object@args[[2]])[1]
  c(lh_length + rh_length - 1, 1)
})

setMethod("sign_from_args", "Conv", function(object) { mul_sign(object@args[[1]], object@args[[2]]) })
setMethod("is_incr", "Conv", function(object, idx) { is_positive(object@args[[1]]) })
setMethod("is_decr", "Conv", function(object, idx) { is_negative(object@args[[1]]) })

Conv.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(conv(arg_objs[[1]], arg_objs[[2]], size), list())
}

setMethod("graph_implementation", "Conv", function(object, arg_objs, size, data = NA_real_) {
  Conv.graph_implementation(arg_objs, size, data)
})

.CumSum <- setClass("CumSum", contains = c("AxisAtom", "AffAtom"))
CumSum <- function(expr, axis = 1) { .CumSum(expr = expr, axis = axis) }
cumsum.Expression <- function(x) { CumSum(expr = Vec(x)) }   # Flatten matrix in column-major order to match R's behavior

setMethod("initialize", "CumSum", function(.Object, ..., expr, axis = 1) {
  .Object@expr <- expr
  .Object@axis <- axis
  callNextMethod(.Object, ..., args = list(.Object@expr))
})

setMethod("to_numeric", "CumSum", function(object, values) { apply(values[[1]], axis, cumsum) })
setMethod("size_from_args", "CumSum", function(object) { size(object@args[[1]]) })

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

.grad.CumSum <- function(object, values) {
  # TODO: This is inefficient
  dim <- dim(values[[1]])[object@axis]
  mat <- matrix(0, nrow = dim, ncol = dim)
  for(i in 1:dim) {
    for(j in 1:(i+1))
      mat[i,j] <- 1
  }
  size <- size(object@args[[1]])
  var <- Variable(size[1], size[2])
  if(object@axis == 2)
    grad <- .grad(MulExpression(mat, var), values)[2]
  else
    grad <- .grad(RMulExpression(var, t(mat)), values)[1]
  list(grad)
}

CumSum.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  # Implicit 0(n) definition:
  # X = Y[:1,:] - Y[1:,:]
  Y <- create_var(size)
  axis <- data[[1]]
  dim <- size[axis]
  diff_mat <- get_diff_mat(dim, axis)
  diff_matt <- create_const(diff_mat, c(dim, dim))
  
  if(axis == 2)
    diff <- mul_expr(diff_mat, Y, size)
  else
    diff <- rmul_expr(Y, diff_mat, size)
  list(Y, list(create_eq(arg_objs[[1]], diff)))
}

setMethod("graph_implementation", "CumSum", function(object, arg_objs, size, data = NA_real_) {
  CumSum.graph_implementation(arg_objs, size, data)
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
  callNextMethod(.Object, ..., args = list(.Object@expr))
})

setMethod("to_numeric", "DiagVec", function(object, values) { diag(as.vector(values[[1]])) })

setMethod("size_from_args", "DiagVec", function(object) {
  rows <- size(object@args[[1]])[1]
  c(rows, rows)
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
  callNextMethod(.Object, ..., args = list(.Object@expr))
})

setMethod("to_numeric", "DiagMat", function(object, values) { diag(values[[1]]) })

setMethod("size_from_args", "DiagMat", function(object) {
  rows <- size(object@args[[1]])[1]
  c(rows, 1)
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
#' @param expr An \S4class{Expression} or numeric constant representing a vector or diagonal matrix.
#' @aliases Diag
#' @export
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

Diff <- function(x, k = 1, axis = 1) {
  x <- as.Constant(x)
  if(axis == 1)
    x <- t(x)
  
  size <- size(x)
  m <- size[1]
  n <- size[2]
  if(k < 0 || k >= m || n != 1)
    stop("Must have k >= 0 and x must be a 1-D vector with < k elements")
  
  d <- x
  for(i in 1:k)
    d <- d[2:length(d)] - d[1:(length(d)-1)]
  
  if(axis == 1)
    t(d)
  else
    d
}

.HStack <- setClass("HStack", contains = "AffAtom")
HStack <- function(...) { .HStack(args = list(...)) }

setMethod("validate_args", "HStack", function(object) {
  arg_cols <- sapply(object@args, function(arg) { size(arg)[1] })
  if(max(arg_cols) != min(arg_cols))
    stop("All arguments to HStack must have the same number of rows")
})

setMethod("to_numeric", "HStack", function(object, values) { Reduce("cbind", values) })
setMethod("size_from_args", "HStack", function(object) {
  arg_cols <- sapply(object@args, function(arg) { size(arg)[2] })
  cols <- sum(arg_cols)
  rows <- size(object@args[[1]])[1]
  c(rows, cols)
})

HStack.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(hstack(arg_objs, size), list())
}

setMethod("graph_implementation", "HStack", function(object, arg_objs, size, data = NA_real_) {
  HStack.graph_implementation(arg_objs, size, data)
})

.Index <- setClass("Index", representation(expr = "Expression", key = "list"), contains = "AffAtom")
Index <- function(expr, key) { .Index(expr = expr, key = key) }

setMethod("initialize", "Index", function(.Object, ..., expr, key) {
  .Object@key <- ku_validate_key(key, size(expr))   # TODO: Double check key validation
  .Object@expr <- expr
  callNextMethod(.Object, ..., args = list(.Object@expr))
})

setMethod("to_numeric", "Index", function(object) {
  ku_slice_mat(values[[1]], object@key$row, object@key$col)
})

setMethod("size_from_args", "Index", function(object) {
  ku_size(object@key, size(object@args[[1]]))
})

setMethod("get_data", "Index", function(object) { list(object@key) })

Index.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  obj <- index(arg_objs[[1]], size, data[[1]])
  list(obj, list())
}

setMethod("graph_implementation", "Index", function(object, arg_objs, size, data = NA_real_) {
  Index.graph_implementation(arg_objs, size, data)
})

Index.get_special_slice <- function(expr, row, col) {
  expr <- as.Constant(expr)
  
  # Order the entries of expr and select them using key.
  expr_size <- size(expr)
  expr_prod <- prod(expr_size)
  
  idx_mat <- seq(expr_prod)
  idx_mat <- matrix(idx_mat, nrow = expr_size[1], ncol = expr_size[2])
  select_mat <- idx_mat[row, col]
  
  if(!is.null(dim(select_mat)))
    final_size <- dim(select_mat)
  else   # Always cast 1-D arrays as column vectors
    final_size <- c(length(select_mat), 1)
  
  # Select the chosen entries from expr.
  select_vec <- as.vector(select_mat)
  identity <- sparseMatrix(i = 1:expr_prod, j = 1:expr_prod, x = rep(1, expr_prod))
  Reshape(identity[select_vec,] %*% Vec(expr), final_size[1], final_size[2])
}

Index.get_index <- function(matrix, constraints, row, col) {
  key <- Key(row, col)
  graph <- Index.graph_implementation(list(matrix), c(1, 1), list(key))
  idx <- graph[[1]]
  idx_constr <- graph[[2]]
  constraints <- c(constraints, idx_constr)
  list(idx = idx, constraints = constraints)
}

Index.block_eq <- function(matrix, block, constraints, row_start, row_end, col_start, col_end) {
  key <- Key(row_start:row_end, col_start:col_end)
  rows <- row_end - row_start
  cols <- col_end - col_start
  if(!all(size(block) == c(rows, cols)))
    stop("Block must have rows = ", rows, " and cols = ", cols)
  graph <- Index.graph_implementation(list(matrix), c(rows, cols), list(key))
  slc <- graph[[1]]
  idx_constr <- graph[[2]]
  constraints <- c(constraints, create_eq(slc, block), idx_constr)
  constraints
}

.Kron <- setClass("Kron", representation(lh_exp = "ConstValORExpr", rh_exp = "ConstValORExpr"), contains = "AffAtom")
Kron <- function(lh_exp, rh_exp) { .Kron(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("initialize", "Kron", function(.Object, ..., lh_exp, rh_exp) {
  .Object@lh_exp <- lh_exp
  .Object@rh_exp <- rh_exp
  callNextMethod(.Object, ..., args = list(.Object@lh_exp, .Object@rh_exp))
})

setMethod("validate_args", "Kron", function(object) {
  if(!is_constant(object@args[[1]]))
    stop("The first argument to Kron must be constant")
})

setMethod("to_numeric", "Kron", function(object, values) {
  kronecker(values[[1]], values[[2]])
})

setMethod("size_from_args", "Kron", function(object) {
  rows <- size(object@args[[1]])[1] * size(object@args[[2]])[1]
  cols <- size(object@args[[1]])[2] * size(object@args[[2]])[2]
  c(rows, cols)
})

setMethod("sign_from_args", "Kron", function(object) { mul_sign(object@args[[1]], object@args[[2]]) })
setMethod("is_incr", "Kron", function(object, idx) { is_positive(object@args[[1]]) })
setMethod("is_decr", "Kron", function(object, idx) { is_negative(object@args[[2]]) })

Kron.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(kron(arg_objs[[1]], arg_objs[[2]], size), list())
}

setMethod("graph_implementation", "Kron", function(object, arg_objs, size, data = NA_real_) {
  Kron.graph_implementation(arg_objs, size, data)
})

.MulElemwise <- setClass("MulElemwise", representation(lh_const = "ConstValORExpr", rh_exp = "ConstValORExpr"), contains = "AffAtom")
MulElemwise <- function(lh_const, rh_exp) { .MulElemwise(lh_const = lh_const, rh_exp = rh_exp) }

setMethod("initialize", "MulElemwise", function(.Object, ..., lh_const, rh_exp) {
  .Object@lh_const <- lh_const
  .Object@rh_exp <- rh_exp
  callNextMethod(.Object, ..., args = list(.Object@lh_const, .Object@rh_exp))
})

setMethod("validate_args", "MulElemwise", function(object) {
  if(!is_constant(object@args[[1]]))
    stop("The first argument to MulElemwise must be constant.")
})

setMethod("to_numeric", "MulElemwise", function(object, values) { values[[1]] * values[[2]] })
setMethod("size_from_args", "MulElemwise", function(object) { 
  sum_shapes(lapply(object@args, function(arg) { size(arg) }))
})

setMethod("sign_from_args", "MulElemwise", function(object) {
  mul_sign(object@args[[1]], object@args[[2]])
})

setMethod("is_incr", "MulElemwise", function(object, idx) { is_positive(object@args[[1]]) })
setMethod("is_decr", "MulElemwise", function(object, idx) { is_negative(object@args[[1]]) })
setMethod("is_quadratic", "MulElemwise", function(object) { is_quadratic(object@args[[2]]) })

MulElemwise.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  if(any(arg_objs[[1]]$size != arg_objs[[2]]$size))
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

setMethod("initialize", "Reshape", function(.Object, ..., expr, rows, cols) {
  .Object@rows <- rows
  .Object@cols <- cols
  .Object@expr <- expr
  callNextMethod(.Object, ..., args = list(.Object@expr))
})

setMethod("validate_args", "Reshape", function(object) {
  old_len <- prod(size(object@args[[1]]))
  new_len <- object@rows * object@cols
  if(old_len != new_len)
    stop(sprintf("Invalid reshape dimensions (%i, %i)", object@rows, object@cols))
})

setMethod("to_numeric", "Reshape", function(object, values) {
  dim(values[[1]]) <- c(object@rows, object@cols)
  values[[1]]
})

setMethod("size_from_args", "Reshape", function(object) { c(object@rows, object@cols) })
setMethod("get_data", "Reshape", function(object) { list(object@rows, object@cols) })

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
#' @slot axis The axis to sum the entries along.
#' @aliases SumEntries
#' @export
.SumEntries <- setClass("SumEntries", contains = c("AxisAtom", "AffAtom"))
SumEntries <- function(expr, axis = NA_real_) { .SumEntries(expr = expr, axis = axis) }

setMethod("to_numeric", "SumEntries", function(object, values) {
  if(is.na(object@axis))
    sum(values[[1]])
  else
    apply(values[[1]], object@axis, sum)
})

SumEntries.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  axis <- data[[1]]
  if(is.na(axis))
    obj <- sum_entries(arg_objs[[1]])
  else if(axis == 2) {
    const_size <- c(arg_objs[[1]]$size[2], 1)
    ones <- create_const(matrix(1, nrow = const_size[1], ncol = const_size[2]), const_size)
    obj <- rmul_expr(arg_objs[[1]], ones, size)
  } else {   # axis == 1
    const_size <- c(1, arg_objs[[1]]$size[1])
    ones <- create_const(matrix(1, nrow = const_size[1], ncol = const_size[2]), const_size)
    obj <- mul_expr(ones, arg_objs[[1]], size)
  }
  list(obj, list())
}

setMethod("graph_implementation", "SumEntries", function(object, arg_objs, size, data = NA_real_) {
  SumEntries.graph_implementation(arg_objs, size, data)
})

sum.Expression <- function(..., na.rm = FALSE) {
  if(na.rm)
    warning("na.rm is unimplemented for Expression objects")
  
  vals <- list(...)
  is_expr <- sapply(vals, function(v) { is(v, "Expression") })
  sum_expr <- lapply(vals[is_expr], function(expr) { SumEntries(expr = expr) })
  if(all(is_expr))
    Reduce("+", sum_expr)
  else {
    sum_num <- sum(sapply(vals[!is_expr], function(v) { sum(v, na.rm = na.rm) }))
    Reduce("+", sum_expr) + sum_num
  }
}

mean.Expression <- function(x, trim = 0, na.rm = FALSE, ...) {
  if(na.rm)
    stop("na.rm is unimplemented for Expression objects")
  if(trim != 0)
    stop("trim is unimplemented for Expression objects")
  SumEntries(expr = x) / prod(size(x))
}

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

setMethod("initialize", "Trace", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., args = list(.Object@expr))
})

setMethod("validate_args", "Trace", function(object) {
  size <- size(object@args[[1]])
  if(size[1] != size[2])
    stop("Argument to trace must be a square matrix")
})

setMethod("to_numeric", "Trace", function(object, values) { sum(diag(values[[1]])) })
setMethod("size_from_args", "Trace", function(object){ c(1, 1) })

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

setMethod("to_numeric", "Transpose", function(object, values) { t(values[[1]]) })

setMethod("size_from_args", "Transpose", function(object) {
  size <- size(object@args[[1]])
  c(size[2], size[1])
})

Transpose.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(transpose(arg_objs[[1]]), list())  
}

setMethod("graph_implementation", "Transpose", function(object, arg_objs, size, data = NA_real_) {
  Transpose.graph_implementation(arg_objs, size, data)
})

.UpperTri <- setClass("UpperTri", representation(expr = "Expression"), contains = "AffAtom")
UpperTri <- function(expr) { .UpperTri(expr = expr) }

setMethod("initialize", "UpperTri", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., args = list(.Object@expr))
})

setMethod("validate_args", "UpperTri", function(object) {
  size <- size(object@args[[1]])
  if(size[1] != size[2])
    stop("Argument to UpperTri must be a square matrix.")
})

setMethod("to_numeric", "UpperTri", function(object, values) {
  # Vectorize the upper triagonal entries
  tridx <- upper.tri(values[[1]], diag = FALSE)
  values[[1]][tridx]
})

setMethod("size_from_args", "UpperTri", function(object) {
  size <- size(object@args[[1]])
  rows <- size[1]
  cols <- size[2]
  c(rows*(cols-1)/2, 1)
})

UpperTri.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(upper_tri(arg_objs[[1]]), list())
}

setMethod("graph_implementation", "UpperTri", function(object, arg_objs, size, data = NA_real_) {
  UpperTri.graph_implementation(arg_objs, size, data)
})

# Flattens the matrix X into a vector in column-major order
Vec <- function(X) {
  X <- as.Constant(X)
  Reshape(expr = X, rows = prod(size(X)), cols = 1)
}

.VStack <- setClass("VStack", contains = "AffAtom")
VStack <- function(...) { .VStack(args = list(...)) }

setMethod("validate_args", "VStack", function(object) {
  arg_cols <- sapply(object@args, function(arg) { size(arg)[2] })
  if(max(arg_cols) != min(arg_cols))
    stop("All arguments to VStack must have the same number of columns")
})

setMethod("to_numeric", "VStack", function(object, values) { Reduce("rbind", values) })
setMethod("size_from_args", "VStack", function(object) {
  cols <- size(object@args[[1]])[2]
  arg_rows <- sapply(object@args, function(arg) { size(arg)[1] })
  rows <- sum(arg_rows)
  c(rows, cols)
})

VStack.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(vstack(arg_objs, size), list())  
}

setMethod("graph_implementation", "VStack", function(object, arg_objs, size, data = NA_real_) {
  VStack.graph_implementation(arg_objs, size, data)
})

Bmat <- function(block_lists) {
  row_blocks <- lapply(block_lists, function(blocks) { .HStack(args = blocks) })
  .VStack(args = row_blocks)
}
