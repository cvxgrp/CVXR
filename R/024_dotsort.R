## CVXPY SOURCE: cvxpy/atoms/dotsort.py

#'
#' The DotSort class.
#'
#' This class represents the value
#'
#' \deqn{\langle sort\left(vec(X)\right), sort\left(vec(W)\right) \rangle},
#'
#' where \eqn{X} is an expression and \eqn{W} is a constant.
#'
#' Both arguments are flattened, i.e., we define \eqn{x = vec(X)} and \eqn{w = vec(W)}.
#' If the length of \eqn{w} is less than the length of \eqn{x}, it is conceptually padded with zeroes.
#' When the length of \eqn{w} is larger than the length of \eqn{x}, an exception is raised.
#'
#' DotSort is a generalization of SumLargest and SumSmallest:
#' SumLargest(X, 3) is equivalent to DotSort(X, c(1,1,1))
#' SumLargest(X, 3.5) is equivalent to DotSort(X, c(1,1,1,0.5))
#' SumSmallest(X, 3) is equivalent to -DotSort(X, c(-1,-1,-1))
#'
#' When the constant argument is not a boolean vector, DotSort can be considered as a weighted sum
#' of \eqn{x}, where the largest weight is assigned to the largest entry in \eqn{x}, etc.
#'
.DotSort <- setClass("DotSort", representation(X = "ConstValORExpr", W = "ConstVal"), contains = "Atom")

#' @param X An \linkS4class{Expression}.
#' @param W A numeric matrix.
#' @rdname DotSort-class
DotSort <- function(X, W) { .DotSort(X = X, W = W) }

setMethod("initialize", "DotSort", function(.Object, ..., X, W) {
  .Object@X <- X
  .Object@W <- W
  callNextMethod(.Object, ..., atom_args = list(.Object@X, .Object@W))
})

#' @describeIn DotSort Check that the arguments are valid.
setMethod("validate_args", "DotSort", function(object) {
  if(!is_constant(object@args[[2]]))
    stop("W must be constant")
  if(size(object@args[[1]]) < size(object@args[[2]]))
    stop("The size of W must be less than or equal to the size of X")
  callNextMethod()
})

#' @param object A \linkS4class{DotSort} object.
#' @param values A list of arguments to the atom.
#' @describeIn DotSort The inner product of the sorted values of vec(X) and the sorted (and potentially padded) values of vec(W).
setMethod("to_numeric", "DotSort", function(object, values) {
  args <- DotSort.get_args_from_values(values)
  x <- args[[1]]
  w_padded <- args[[2]]
  return(as.numeric(sort(x) %*% sort(w_padded)))
})

#' @param values A list of numeric values for the arguments
#' @describeIn DotSort Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "DotSort", function(object, values) {
  args <- DotSort.get_args_from_values(values)
  x <- args[[1]]
  w_padded <- args[[2]]
  n <- length(x)
  sorted_w <- sort(w_padded)
  return(list(sparseMatrix(i = indices, j = rep(1, n), x = sorted_w, dims = c(n, 1))))
})

#' @describeIn DotSort The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "DotSort", function(object) { c(1,1) })

#' @describeIn DotSort The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "DotSort", function(object) {
  # Same as argument.
  x_pos <- is_nonneg(object@args[[1]])
  x_neg <- is_nonpos(object@args[[1]])

  w_pos <- is_nonneg(object@args[[2]])
  w_neg <- is_nonpos(object@args[[2]])

  is_positive <- (x_pos && w_pos) || (x_neg && w_neg)
  is_negative <- (x_neg && w_pos) || (x_pos && w_neg)

  return(c(is_positive, is_negative))
})

#' @describeIn DotSort Is the atom convex?
setMethod("is_atom_convex", "DotSort", function(object) {
  if(dpp_scope_active()) {   # TODO: Figure out how to save DPP status as global parameter.
    # DotSort is convex under DPP if W is parameter affine.
    X <- object@args[[1]]
    W <- object@args[[2]]
    return(is_constant(X) || is_param_affine(W))
  } else
    return(TRUE)
})

#' @describeIn DotSort Is the atom concave?
setMethod("is_atom_concave", "DotSort", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn DotSort Is the atom weakly increasing in the index?
setMethod("is_incr", "DotSort", function(object, idx) { is_nonneg(object@args[[2]]) })

#' @describeIn DotSort Is the atom weakly decreasing in the index?
setMethod("is_decr", "DotSort", function(object, idx) { is_nonpos(object@args[[2]]) })

#' @describeIn DotSort Empty list because W is stored as an argument.
setMethod("get_data", "DotSort", function(object) { list() })

DotSort.get_args_from_values <- function(values) {
  x <- as.numeric(t(values[[1]]))
  w <- as.numeric(t(values[[2]]))

  w_padded <- rep(0, length(x))
  w_padded[1:length(w)] <- w
  return(list(x, w_padded))
}
