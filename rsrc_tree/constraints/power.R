## CVXPY SOURCE: cvxpy/expression/constraints/power.py
#'
#' The PowCone3D class.
#'
#' This class represents a collection of 3D power cone constraints
#'
#' \deqn{x[i]^alpha[i] * y[i]^(1-alpha[i]) \geq |z[i]| \text{for all} i, x \geq 0, y \geq 0}
#'
#' If the parameter alpha is a scalar, it will be promoted to a vector matching
#' the (common) sizes of x, y, z. The numeric value of alpha (or its components,
#' in the vector case) must be a number in the open interval (0, 1).
#'
#' We store flattened representations of the arguments (x, y, z, and alpha) as
#' Expression objects. We construct dual variables with respect to these
#' flattened representations.
#'
#' @slot x An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{x}.
#' @slot y An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{y}.
#' @slot z An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{z}.
#' @slot alpha An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{\alpha}. Must be in the open interval (0, 1).
#' @name PowCone3D-class
#' @aliases PowCone3D
#' @rdname PowCone3D-class
.PowCone3D <- setClass("PowCone3D", representation(x = "ConstValORExpr", y = "ConstValORExpr", z = "ConstValORExpr", alpha = "ConstValORExpr"), contains = "Constraint")

#' @param x An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{x}.
#' @param y An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{y}.
#' @param z An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{z}.
#' @param alpha An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{\alpha}. Must be in the open interval (0, 1).
#' @rdname PowCone3D-class
PowCone3D <- function(x, y, z, alpha, constr_id = NA_integer_) { .PowCone3D(x = x, y = y, z = z, alpha = alpha, constr_id = constr_id) }

setMethod("initialize", "PowCone3D", function(.Object, ..., x, y, z, alpha) {
  .Object@x <- as.Constant(x)
  .Object@y <- as.Constant(y)
  .Object@z <- as.Constant(z)
  for(val in list(.Object@x, .Object@y, .Object@z)) {
    if(!(is_affine(val) && is_real(val)))
      stop("All arguments must be affine and real")
  }

  alpha <- as.Constant(alpha)
  if(is_scalar(alpha))
    alpha <- Promote(alpha, dim(.Object@x))
  .Object@alpha <- alpha

  alpha_val <- value(.Object@alpha)
  if(any(alpha_val <= 0) || any(alpha_val >= 1))
    stop("alpha must have entries in the open interval (0, 1)")
  arg_dims <- list(dim(.Object@x), dim(.Object@y), dim(.Object@z), dim(.Object@alpha))
  for(i in 2:length(arg_dims)) {
    s <- arg_dims[[i]]
    if(any(arg_dims[[1]] != s))
      stop("All arguments must have the same dimensions")
  }
  callNextMethod(.Object, ..., args = list(.Object@x, .Object@y, .Object@z))
})

#' @rdname PowCone3D-class
setMethod("as.character", "PowCone3D", function(x) {
  paste("PowCone3D(", as.character(x@x), ", ", as.character(x@y), ", ", as.character(x@z), ", ", as.character(x@alpha), ")", sep = "")
})

#' @describeIn PowCone3D A \linkS4class{Expression} representing the residual of the constraint.
setMethod("residual", "PowCone3D", function(object) {
  # TODO: The projection should be implemented directly.
  if(is.na(value(object@x)) || is.na(value(object@y)) || is.na(value(object@z)))
    return(NA_real_)

  x <- new("Variable", dim = dim(object@x))
  y <- new("Variable", dim = dim(object@y))
  z <- new("Variable", dim = dim(object@z))
  constr <- list(PowCone3D(x, y, z, object@alpha))
  obj <- Minimize(Norm2(HStack(x, y, z) - HStack(value(object@x), value(object@y), value(object@z))))
  problem <- Problem(obj, constr)
  result <- solve(problem, solver = "SCS", eps = 1e-8)
  return(result$value)
})

#' @describeIn PowCone3D Information needed to reconstruct the object aside from the args.
setMethod("get_data", "PowCone3D", function(object) { list(object@alpha, id(object)) })

#' @describeIn PowCone3D A logical value indicating whether the constraint is imaginary.
setMethod("is_imag", "PowCone3D", function(object) { FALSE })

#' @describeIn PowCone3D A logical value indicating whether the constraint is complex.
setMethod("is_complex", "PowCone3D", function(object) { FALSE })

#' @describeIn PowCone3D The number of elementwise cones.
setMethod("num_cones", "PowCone3D", function(object) { size(object@x) })

#' @describeIn PowCone3D The dimensions of the second-order cones.
setMethod("cone_sizes", "PowCone3D", function(object) { rep(3, num_cones(object)) })

#' @describeIn PowCone3D The constraint is DCP if the constrained expression is affine.
setMethod("is_dcp", "PowCone3D", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping
    args_ok <- all(sapply(object@args, is_affine))
    exps_ok <- !is(object@alpha, "Parameter")
    return(args_ok && exps_ok)
  }
  return(all(sapply(object@args, is_affine)))
})

#' @describeIn PowCone3D Is the constraint DGP?
setMethod("is_dgp", "PowCone3D", function(object, dpp = FALSE) { FALSE })

#' @describeIn PowCone3D Is the constraint DQCP?
setMethod("is_dqcp", "PowCone3D", function(object) { is_dcp(object) })

#' @describeIn PowCone3D The dimensions of the constrained expression.
setMethod("dim", "PowCone3D", function(x) {
  s <- c(3, dim(x@x))
  # Note: this can be a 3-tuple of ndim(x) == 2.
  return(s)
})

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn PowCone3D Replaces the dual values of a second-order cone constraint.
setReplaceMethod("dual_value", "PowCone3D", function(object, value) {
  # TODO: figure out why the reshaping has to be done differently,
  # compared to ExpCone constraints.
  value <- t(matrix(t(value), ncol = 3, byrow = FALSE))
  dv0 <- matrix(value[1,], nrow = nrow(object@x), ncol = ncol(object@x))
  dv1 <- matrix(value[2,], nrow = nrow(object@y), ncol = ncol(object@y))
  dv2 <- matrix(value[3,], nrow = nrow(object@z), ncol = ncol(object@z))

  value(object@dual_variables[[1]]) <- dv0
  value(object@dual_variables[[2]]) <- dv1
  value(object@dual_variables[[3]]) <- dv2
  return(object)
})

#'
#' The PowConeND class.
#'
#' This class represents a collection of N-dimensional power cone constraints
#' that is mathematically equivalent to the following code snippet:
#'
#' \code{apply(W^alpha, axis, prod) >= abs(z)},
#' W >= 0
#'
#' All arguments must be Expression-like, and z must satisfy ndim(z) <= 1. The
#' rows (resp. columns) of alpha must sum to 1 when axis = 1 (resp. axis = 2).
#'
#' Note: unlike PowCone3D, we make no attempt to promote alpha to the
#' appropriate shape. The dimensions of W and alpha must match exactly.
#'
#' Note: Dual variables are not currently implemented for this type of constraint.
#'
#' @slot W An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{W}.
#' @slot z An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{z}.
#' @slot alpha An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{\alpha}. Must be in the open interval (0, 1).
#' @slot axis The dimension along which to constrain: \code{1} indicates rows, and \code{2} indicates columns. The default is \code{2}.
#' @name PowConeND-class
#' @aliases PowConeND
#' @rdname PowConeND-class
.PowConeND <- setClass("PowConeND", representation(W = "ConstValORExpr", z = "ConstValORExpr", alpha = "ConstValORExpr", axis = "numeric"), prototype(axis = 2),
                       validity = function(object) {
                                      if(length(axis) > 1 || (axis != 1 && axis != 2))
                                        stop("[PowConeND: axis] axis must be either 1 (rows) or 2 (columns)")
                                      return(TRUE)
                                    }, contains = "Constraint")

#' @param W An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{W}.
#' @param z An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{z}.
#' @param alpha An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{\alpha}. Must be in the open interval (0, 1).
#' @param axis The dimension along which to constrain: \code{1} indicates rows, and \code{2} indicates columns. The default is \code{2}.
#' @rdname PowConeND-class
PowConeND <- function(W, z, alpha, axis = 2, constr_id = NA_integer_) { .PowConeND(W = W, z = z, alpha = alpha, axis = axis, constr_id = constr_id) }

setMethod("initialize", "PowConeND", function(.Object, ..., W, z, alpha, axis = 2) {
  W <- as.Constant(W)
  if(!(is_real(W) && is_affine(W)))
    stop("Invalid first argument; W must be affine and real.")

  z <- as.Constant(z)
  # if(!(ndim(z) <= 1 || (ndim(z) == 2 && ncol(z) == 1)) || !(is_real(z) && is_affine(z)))
  if(ndim(z) > 1 || !(is_real(z) && is_affine(z)))
    stop("Invalid second argument. z must be affine, real, and have at most one ndim(z) <= 1.")

  # Check z has one entry per cone.
  if((ndim(W) <= 1 && size(z) > 1) ||
     (ndim(W) == 2 && size(z) != dim(W)[axis]) ||
     (ndim(W) == 1 && axis == 1))
    stop("Argument dimensions and axis are incompatible")

  axis_opp <- ifelse(axis == 1, 2, 1)
  if(ndim(W) == 2 && dim(W)[axis_opp] <= 1)
    stop("PowConeND requires left-hand-side to have at least two terms.")

  alpha <- as.Constant(alpha)
  if(any(dim(alpha) != dim(W)))
    stop("W and alpha dimensions must be equal")
  if(any(value(alpha) <= 0))
    stop("Argument alpha must be entry-wise positive.")
  if(any(abs(1 - apply(value(alpha), axis_opp, sum)) > 1e-6))
    stop("Argument alpha must sum to 1 along specified axis.")

  .Object@W <- W
  .Object@z <- z
  .Object@alpha <- alpha
  .Object@axis <- axis
  if(ndim(z) == 0)
    z <- flatten(z)

  callNextMethod(.Object, ..., args = list(W, z))
})

#' @rdname PowConeND-class
setMethod("as.character", "PowConeND", function(x) {
  paste("PowConeND(", as.character(x@x), ", ", as.character(x@W), ", ", as.character(x@z), ", ", as.character(x@alpha), ")", sep = "")
})

#' @describeIn PowConeND A logical value indicating whether the constraint is imaginary.
setMethod("is_imag", "PowConeND", function(object) { FALSE })

#' @describeIn PowConeND A logical value indicating whether the constraint is complex.
setMethod("is_complex", "PowConeND", function(object) { FALSE })

#' @describeIn PowConeND Information needed to reconstruct the object aside from the args.
setMethod("get_data", "PowConeND", function(object) { list(object@alpha, object@axis, id(object)) })

#' @describeIn PowConeND A \linkS4class{Expression} representing the residual of the constraint.
setMethod("residual", "PowConeND", function(object) {
  # TODO: The projection should be implemented directly.
  if(is.na(value(object@W)) || is.na(value(object@z)))
    return(NA_real_)

  W <- new("Variable", dim = dim(object@W))
  z <- new("Variable", dim = dim(object@z))
  constr <- list(PowConeND(W, z, object@alpha, axis = object@axis))
  obj <- Minimize(Norm2(HStack(flatten(W), flatten(z)) -
                        HStack(value(flatten(object@W)), value(flatten(object@z)))))
  problem <- Problem(obj, constr)
  result <- solve(problem, solver = "SCS", eps = 1e-8)
  return(result$value)
})

#' @describeIn PowConeND The number of elementwise cones.
setMethod("num_cones", "PowConeND", function(object) { size(object@z) })

#' @describeIn PowConeND The number of entries in the combined cones.
setMethod("size", "PowConeND", function(object) {
  axis_opp <- ifelse(object@axis == 1, 2, 1)
  cone_size <- 1 + dim(object@args[[1]])[axis_opp]
  return(cone_size * num_cones(object))
})

#' @describeIn PowConeND The dimensions of the second-order cones.
setMethod("cone_sizes", "PowConeND", function(object) {
  axis_opp <- ifelse(object@axis == 1, 2, 1)
  cone_size <- 1 + dim(object@args[[1]])[axis_opp]
  rep(cone_size, num_cones(object))
})

#' @describeIn PowConeND The constraint is DCP if the constrained expression is affine.
setMethod("is_dcp", "PowConeND", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping
    args_ok <- (is_affine(object@args[[1]]) && is_affine(object@args[[2]]))
    exps_ok <- !is(object@alpha, "Parameter")
    return(args_ok && exps_ok)
  }
  return(TRUE)
})

#' @describeIn PowConeND Is the constraint DGP?
setMethod("is_dgp", "PowConeND", function(object, dpp = FALSE) { FALSE })

#' @describeIn PowConeND Is the constraint DQCP?
setMethod("is_dqcp", "PowConeND", function(object) { is_dcp(object) })

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn PowConeND Replaces the dual values of a second-order cone constraint.
setReplaceMethod("dual_value", "PowConeND", function(object, value) {
  stop("Unimplemented")
})
