## CVXPY SOURCE: cvxpy/atoms/affine/index.py

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
  vec_arg <- lu.reshape(arg, c(arg_size, 1))
  mul_mat <- id_mat[select_vec]
  mul_const <- lu.create_const(mul_mat, dim(mul_mat), sparse = TRUE)
  mul_expr <- lu.mul_expr(mul_const, vec_arg, c(nrow(mul_mat), 1))
  obj <- lu.reshape(mul_expr, final_dim)
  return(list(obj, list()))
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Index The graph implementation of the atom.
setMethod("graph_implementation", "SpecialIndex", function(object, arg_objs, dim, data = NA_real_) {
  SpecialIndex.graph_implementation(arg_objs, dim, data)
})
