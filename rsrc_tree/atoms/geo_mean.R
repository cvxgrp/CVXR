## CVXPY SOURCE: cvxpy/atoms/geo_mean.py

#'
#' The GeoMean class.
#'
#' This class represents the (weighted) geometric mean of vector \eqn{x} with optional powers given by \eqn{p}.
#'
#' \deqn{\left(x_1^{p_1} \cdots x_n^{p_n} \right)^{\frac{1}{\mathbf{1}^Tp}}}
#'
#' The geometric mean includes an implicit constraint that \eqn{x_i \geq 0} whenever \eqn{p_i > 0}. If \eqn{p_i = 0, x_i} will be unconstrained.
#' The only exception to this rule occurs when \eqn{p} has exactly one nonzero element, say \eqn{p_i}, in which case \code{GeoMean(x,p)} is equivalent to \eqn{x_i} (without the nonnegativity constraint).
#' A specific case of this is when \eqn{x \in \mathbf{R}^1}.
#'
#' @slot x An \linkS4class{Expression} or numeric vector.
#' @slot p (Optional) A vector of weights for the weighted geometric mean. The default is a vector of ones, giving the \strong{unweighted} geometric mean \eqn{x_1^{1/n} \cdots x_n^{1/n}}.
#' @slot max_denom (Optional) The maximum denominator to use in approximating \code{p/sum(p)} with \code{w}. If \code{w} is not an exact representation, increasing \code{max_denom} may offer a more accurate representation, at the cost of requiring more convex inequalities to represent the geometric mean. Defaults to 1024.
#' @slot w (Internal) A list of \code{bigq} objects that represent a rational approximation of \code{p/sum(p)}.
#' @slot approx_error (Internal) The error in approximating \code{p/sum(p)} with \code{w}, given by \eqn{\|p/\mathbf{1}^Tp - w\|_{\infty}}.
#' @name GeoMean-class
#' @aliases GeoMean
#' @importClassesFrom gmp bigq bigz
#' @rdname GeoMean-class
.GeoMean <- setClass("GeoMean", representation(x = "ConstValORExpr", p = "numeric", max_denom = "numeric",
                                               w = "bigq", w_dyad = "bigq", approx_error = "numeric", tree = "Rdict",
                                               cone_lb = "numeric", cone_num = "numeric", cone_num_over = "numeric"),
                                prototype(p = NA_real_, max_denom = 1024), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric vector.
#' @param p (Optional) A vector of weights for the weighted geometric mean. The default is a vector of ones, giving the \strong{unweighted} geometric mean \eqn{x_1^{1/n} \cdots x_n^{1/n}}.
#' @param max_denom (Optional) The maximum denominator to use in approximating \code{p/sum(p)} with \code{w}. If \code{w} is not an exact representation, increasing \code{max_denom} may offer a more accurate representation, at the cost of requiring more convex inequalities to represent the geometric mean. Defaults to 1024.
#' @rdname GeoMean-class
GeoMean <- function(x, p = NA_real_, max_denom = 1024) { .GeoMean(x = x, p = p, max_denom  = max_denom) }

setMethod("initialize", "GeoMean", function(.Object, ..., x, p = NA_real_, max_denom = 1024) {
  .Object@x <- x
  .Object@max_denom <- max_denom
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@x), validate = FALSE)

  x <- .Object@args[[1]]
  if(is_vector(x))
    n <- ifelse(ndim(x) == 0, 1, max(dim(x)))
  else
    stop("x must be a row or column vector.")

  if(any(is.na(p)))
    p <- rep(1, n)
  .Object@p <- p

  if(length(.Object@p) != n)
    stop("x and p must have the same number of elements.")

  if(any(.Object@p < 0) || sum(.Object@p) <= 0)
    stop("powers must be nonnegative and not all zero.")

  frac <- fracify(.Object@p, .Object@max_denom)
  .Object@w <- frac[[1]]
  .Object@w_dyad <- frac[[2]]
  .Object@approx_error <- approx_error(.Object@p, .Object@w)

  .Object@tree <- decompose(.Object@w_dyad)

  # known lower bound on number of cones needed to represent w_dyad
  .Object@cone_lb <- lower_bound(.Object@w_dyad)

  # number of cones used past known lower bound
  .Object@cone_num_over <- over_bound(.Object@w_dyad, .Object@tree)

  # number of cones used
  .Object@cone_num <- .Object@cone_lb + .Object@cone_num_over
  validObject(.Object)
  .Object
})

#' @param object A \linkS4class{GeoMean} object.
#' @param values A list of arguments to the atom.
#' @describeIn GeoMean The (weighted) geometric mean of the elements of \code{x}.
setMethod("to_numeric", "GeoMean", function(object, values) {
  values <- as.vector(values[[1]])
  val <- 1.0
  for(idx in 1:length(values)) {
    x <- values[[idx]]
    p <- object@w[idx]
    val <- val * Rmpfr::mpfr(x, Rmpfr::getPrec(x))^p
  }
  return(gmp::asNumeric(val))   # TODO: Handle mpfr objects in the backend later
})

#' @describeIn GeoMean Returns constraints describing the domain of the node
setMethod(".domain", "GeoMean", function(object) { list(object@args[[1]][object@w > 0] >= 0) })

#' @param values A list of numeric values for the arguments
#' @describeIn GeoMean Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "GeoMean", function(object, values) {
  x <- as.matrix(values[[1]])
  # No special case when only one non-zero weight
  w_arr <- as.double(object@w)   # TODO: I'm casting bigq/bigz to double to construct Matrix properly.
  # Outside domain
  if(any(x[w_arr > 0] <= 0))
    return(list(NA_real_))
  else {
    D <- w_arr/as.vector(x) * to_numeric(object, values)
    return(list(Matrix(D, sparse = TRUE)))
  }
})

#' @describeIn GeoMean The name and arguments of the atom.
setMethod("name", "GeoMean", function(x) {
  vals <- paste(sapply(x@w, as.character), collapse = ", ")
  paste("GeoMean(", name(x@args[[1]]), ", (", vals, "))", sep = "")
})

#' @describeIn GeoMean The atom is a scalar.
setMethod("dim_from_args", "GeoMean", function(object) { c(1,1) })

#' @describeIn GeoMean The atom is non-negative.
setMethod("sign_from_args", "GeoMean", function(object) { c(TRUE, FALSE) })

#' @describeIn GeoMean The atom is not convex.
setMethod("is_atom_convex", "GeoMean", function(object) { FALSE })

#' @describeIn GeoMean The atom is concave.
setMethod("is_atom_concave", "GeoMean", function(object) { TRUE })

#' @describeIn GeoMean Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "GeoMean", function(object) { TRUE })

#' @describeIn GeoMean Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "GeoMean", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn GeoMean The atom is weakly increasing in every argument.
setMethod("is_incr", "GeoMean", function(object, idx) { TRUE })

#' @describeIn GeoMean The atom is not weakly decreasing in any argument.
setMethod("is_decr", "GeoMean", function(object, idx) { FALSE })

#' @describeIn GeoMean Returns \code{list(w, dyadic completion, tree of dyads)}.
setMethod("get_data", "GeoMean", function(object) { list(object@w, object@w_dyad, object@tree) })

#' @param args An optional list that contains the arguments to reconstruct the atom. Default is to use current arguments of the atom.
#' @param id_objects Currently unused.
#' @describeIn GeoMean Returns a shallow copy of the GeoMean atom
setMethod("copy", "GeoMean", function(object, args = NULL, id_objects = list()) {
  if(is.null(args))
    args <- object@args

  copy <- do.call(class(object), args)
  data <- get_data(object)
  copy@w <- data[[1]]
  copy@w_dyad <- data[[2]]
  copy@tree <- data[[3]]

  copy@approx_error <- object@approx_error
  copy@cone_lb <- object@cone_lb
  copy@cone_num_over <- object@cone_num_over
  copy@cone_num <- object@cone_num
  copy
})

