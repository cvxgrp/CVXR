## CVXPY SOURCE: cvxpy/expression/consraints/exponential.py
#'
#' The ExpCone class.
#'
#' This class represents a reformulated exponential cone constraint operating elementwise on \eqn{a, b, c}.
#'
#' Original cone:
#' \deqn{
#' K = \{(x,y,z) | y > 0, ye^{x/y} \leq z\} \cup \{(x,y,z) | x \leq 0, y = 0, z \geq 0\}
#' }
#' Reformulated cone:
#' \deqn{
#' K = \{(x,y,z) | y, z > 0, y\log(y) + x \leq y\log(z)\} \cup \{(x,y,z) | x \leq 0, y = 0, z \geq 0\}
#' }
#'
#' @slot x The variable \eqn{x} in the exponential cone.
#' @slot y The variable \eqn{y} in the exponential cone.
#' @slot z The variable \eqn{z} in the exponential cone.
#' @name ExpCone-class
#' @aliases ExpCone
#' @rdname ExpCone-class
.ExpCone <- setClass("ExpCone", representation(x = "ConstValORExpr", y = "ConstValORExpr", z = "ConstValORExpr"), contains = "Constraint")

#' @param x The variable \eqn{x} in the exponential cone.
#' @param y The variable \eqn{y} in the exponential cone.
#' @param z The variable \eqn{z} in the exponential cone.
#' @param constr_id (Optional) A numeric value representing the constraint ID.
#' @rdname ExpCone-class
## #' @export
ExpCone <- function(x, y, z, constr_id = NA_integer_) { .ExpCone(x = x, y = y, z = z, constr_id = constr_id) }

setMethod("initialize", "ExpCone", function(.Object, ..., x, y, z) {
  .Object@x <- as.Constant(x)
  .Object@y <- as.Constant(y)
  .Object@z <- as.Constant(z)
  args <- list(.Object@x, .Object@y, .Object@z)
  for(val in args) {
    if(!(is_affine(val) && is_real(val)))
      stop("All arguments must be affine and real")
  }
  xs <- dim(.Object@x)
  ys <- dim(.Object@y)
  zs <- dim(.Object@z)
  if(!all(xs == ys) || !all(xs == zs))
    stop("All arguments must have the same shapes. Provided arguments have shapes ", xs, ", ", ys, ", and ", zs)
  callNextMethod(.Object, ..., args = list(.Object@x, .Object@y, .Object@z))
})

setMethod("show", "ExpCone", function(object) {
  print(paste("ExpCone(", as.character(object@x), ", ", as.character(object@y), ", ", as.character(object@z), ")", sep = ""))
})

#' @rdname ExpCone-class
setMethod("as.character", "ExpCone", function(x) {
  paste("ExpCone(", as.character(x@x), ", ", as.character(x@y), ", ", as.character(x@z), ")", sep = "")
})

#' @param object A \linkS4class{ExpCone} object.
#' @describeIn ExpCone The size of the \code{x} argument.
setMethod("residual", "ExpCone", function(object) {
  # TODO: The projection should be implemented directly.
  if(any(is.na(value(object@x))) || any(is.na(value(object@y))) || any(is.na(value(object@z))))
    return(NA_real_)
  # x <- Variable(dim(object@x))
  # y <- Variable(dim(object@y))
  # z <- Variable(dim(object@z))
  x <- new("Variable", dim = dim(object@x))
  y <- new("Variable", dim = dim(object@y))
  z <- new("Variable", dim = dim(object@z))
  constr <- list(ExpCone(x, y, z))
  obj <- Minimize(Norm2(HStack(x, y, z) - HStack(value(object@x), value(object@y), value(object@z))))
  prob <- Problem(obj, constr)
  result <- solve(prob)
  return(result$value)
})

#' @describeIn ExpCone The number of entries in the combined cones.
setMethod("size", "ExpCone", function(object) { 3 * num_cones(object) })

#' @describeIn ExpCone The number of elementwise cones.
setMethod("num_cones", "ExpCone", function(object) { size(object@x) })

setMethod("as_quad_approx", "ExpCone", function(object, m, k) {
  if(as.integer(m) != m)
    stop("m must be an integer")
  if(as.integer(k) != k)
    stop("k must be an integer")
  return(RelEntrConeQuad(object@y, object@z, -object@x, m, k))
})

#' @describeIn ExpCone The dimensions of the exponential cones.
setMethod("cone_sizes", "ExpCone", function(object) { rep(3, num_cones(object)) })

#' @describeIn ExpCone An exponential constraint is DCP if each argument is affine.
setMethod("is_dcp", "ExpCone", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping
    return(all(sapply(object@args, is_affine)))
  }
  return(all(sapply(object@args, is_affine)))
})

#' @describeIn ExpCone Is the constraint DGP?
setMethod("is_dgp", "ExpCone", function(object, dpp = FALSE) { FALSE })

#' @describeIn ExpCone Is the constraint DQCP?
setMethod("is_dqcp", "ExpCone", function(object) { is_dcp(object) })

setMethod("dim", "ExpCone", function(x) { c(3, dim(x@x)) })

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn ExpCone Replaces the dual values of an exponential cone constraint.
setReplaceMethod("dual_value", "ExpCone", function(object, value) {
  # TODO: verify that reshaping below works correctly
  value <- t(matrix(t(value), nrow = 3, byrow = FALSE))
  dv0 <- matrix(value[, 1], nrow = nrow(object@x), ncol = ncol(object@x))
  dv1 <- matrix(value[, 2], nrow = nrow(object@y), ncol = ncol(object@y))
  dv2 <- matrix(value[, 3], nrow = nrow(object@z), ncol = ncol(object@z))

  value(object@dual_variables[[1]]) <- dv0
  value(object@dual_variables[[2]]) <- dv1
  value(object@dual_variables[[3]]) <- dv2
  return(object)
})

#'
#' The RelEntrConeQuad class.
#'
#' This class represents an approximate construction of the scalar relative entropy cone.
#'
#' \deqn{
#'  K_{re}=\\text{cl}\\{(x,y,z)\\in\\mathbb{R}_{++} \\times \\mathbb{R}_{++}\\times\\mathbb{R}_{++}\\:x\\log(x/y)\\leq z\\}
#' }
#'
#' Since the above definition is very similar to the ExpCone, we provide a conversion method.
#'
#' More details on the approximation can be found in Theorem-3 on page-10 in the paper:
#' Semidefinite Approximations of the Matrix Logarithm.
#'
#' @slot x The variable \eqn{x} in the (approximate) scalar relative entropy cone.
#' @slot y The variable \eqn{y} in the (approximate) scalar relative entropy cone.
#' @slot z The variable \eqn{z} in the (approximate) scalar relative entropy cone.
#' @slot m An integer directly related to the number of generated nodes for the quadrature approximation used in the algorithm.
#' @slot k An integer controlling the approximation.
#' @name RelEntrConeQuad-class
#' @aliases RelEntrConeQuad
#' @rdname RelEntrConeQuad-class
.RelEntrConeQuad <- setClass("RelEntrConeQuad", representation(x = "ConstValORExpr", y = "ConstValORExpr", z = "ConstValORExpr", m = "numeric", k = "numeric"),
                             validity = function(object) {
                               if(as.integer(object@m) != object@m)
                                 stop("[RelEntrConeQuad: m] The argument m must be an integer")
                               if(as.integer(object@k) != object@k)
                                 stop("[RelEntrConeQuad: k] The argument k must be an integer")
                               return(TRUE)
                              }, contains = "Constraint")

#' @param x The variable \eqn{x} in the (approximate) scalar relative entropy cone.
#' @param y The variable \eqn{y} in the (approximate) scalar relative entropy cone.
#' @param z The variable \eqn{z} in the (approximate) scalar relative entropy cone.
#' @param m An integer directly related to the number of generated nodes for the quadrature approximation used in the algorithm.
#' @param k An integer controlling the approximation.
#' @param constr_id (Optional) A numeric value representing the constraint ID.
#' @rdname RelEntrConeQuad-class
## #' @export
RelEntrConeQuad <- function(x, y, z, m, k, constr_id = NA_integer_) { .RelEntrConeQuad(x = x, y = y, z = z, m = m, k = k, constr_id = constr_id) }

setMethod("initialize", "RelEntrConeQuad", function(.Object, ..., x, y, z, m, k) {
  .Object@x <- as.Constant(x)
  .Object@y <- as.Constant(y)
  .Object@z <- as.Constant(z)
  args <- list(.Object@x, .Object@y, .Object@z)
  for(val in args) {
    if(!(is_affine(val) && is_real(val)))
      stop("All arguments must be affine and real")
  }
  .Object@m <- m
  .Object@k <- k
  xs <- dim(.Object@x)
  ys <- dim(.Object@y)
  zs <- dim(.Object@z)
  if(!all(xs == ys) || !all(xs == zs))
    stop("All arguments must have the same shapes. Provided arguments have shapes ", xs, ", ", ys, ", and ", zs)
  callNextMethod(.Object, ..., args = list(.Object@x, .Object@y, .Object@z))
})

#' @describeIn RelEntrConeQuad Information needed to reconstruct the object aside from the args.
setMethod("get_data", "RelEntrConeQuad", function(object) { list(object@m, object@k, id(object)) })

setMethod("show", "RelEntrConeQuad", function(object) {
  print(paste("RelEntrConeQuad(", as.character(object@x), ", ", as.character(object@y), ", ", as.character(object@z), ", ", object@m, ", ", object@k, ")", sep = ""))
})

#' @rdname RelEntrConeQuad-class
setMethod("as.character", "RelEntrConeQuad", function(x) {
  paste("RelEntrConeQuad(", as.character(x@x), ", ", as.character(x@y), ", ", as.character(x@z), ", ", x@m, ", ", x@k, ")", sep = "")
})

setMethod("residual", "RelEntrConeQuad", function(object) {
  # TODO: The projection should be implemented directly.
  if(is.na(value(object@x)) || is.na(value(object@y)) || is.na(value(object@z)))
    return(NA_real_)

  x <- new("Variable", dim = dim(object@x))
  y <- new("Variable", dim = dim(object@y))
  z <- new("Variable", dim = dim(object@z))
  constr <- list(RelEntrConeQuad(x, y, z, object@m, object@k))
  obj <- Minimize(Norm2(HStack(x, y, z) - HStack(value(object@x), value(object@y), value(object@z))))
  problem <- Problem(obj, constr)
  result <- solve(problem)
  return(result$value)
})

#' @describeIn RelEntrConeQuad The number of entries in the combined cones.
setMethod("size", "RelEntrConeQuad", function(object) { 3*num_cones(object) })

#' @describeIn RelEntrConeQuad The number of elementwise cones.
setMethod("num_cones", "RelEntrConeQuad", function(object) { size(object@x) })

#' @describeIn RelEntrConeQuad The dimensions of the exponential cones.
setMethod("cone_sizes", "RelEntrConeQuad", function(object) { rep(3, num_cones(object)) })

#' @describeIn RelEntrConeQuad An exponential constraint is DCP if each argument is affine.
setMethod("is_dcp", "RelEntrConeQuad", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping
    return(all(sapply(object@args, is_affine)))
  }
  return(all(sapply(object@args, is_affine)))
})

#' @describeIn RelEntrConeQuad Is the constraint DGP?
setMethod("is_dgp", "RelEntrConeQuad", function(object, dpp = FALSE) { FALSE })

#' @describeIn RelEntrConeQuad Is the constraint DQCP?
setMethod("is_dqcp", "RelEntrConeQuad", function(object) { is_dcp(object) })

setMethod("dim", "RelEntrConeQuad", function(x) { c(3, dim(x@x)) })

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn RelEntrConeQuad Replaces the dual values of an exponential cone constraint.
setReplaceMethod("dual_value", "RelEntrConeQuad", function(object, value) {
  stop("Unimplemented")
})

#'
#' The OpRelEntrConeQuad class.
#'
#' This class represents an approximate construction of the scalar relative entropy cone.
#'
#' \deqn{
#'  K_{re}^n=\\text{cl}\\{(X,Y,T)\\in\\mathbb{H}^n_{++} \\times \\mathbb{H}^n_{++}\\times\\mathbb{H}^n_{++}\\:D_{\\text{op}}\\succeq T\\}
#' }
#'
#' This approximation uses \eqn{m + k} semidefinite constraints.
#'
#' More details on the approximation can be found in Theorem-3 on page-10 in the paper:
#' Semidefinite Approximations of the Matrix Logarithm.
#'
#' @slot X The variable \eqn{X} in the (approximate) operator relative entropy cone.
#' @slot Y The variable \eqn{Y} in the (approximate) operator relative entropy cone.
#' @slot Z The variable \eqn{Z} in the (approximate) operator relative entropy cone.
#' @slot m A positive integer that controls the number of quadrature nodes used in a local approximation of the matrix logarithm. Increasing this value results in better local approximations, but does not significantly expand the region of inputs for which the approximation is effective.
#' @slot k A positive integer that sets the number of scaling points about which the quadrature approximation is performed. Increasing this value will expand the region of inputs over which the approximation is effective.
#' @name OpRelEntrConeQuad-class
#' @aliases OpRelEntrConeQuad
#' @rdname OpRelEntrConeQuad-class
.OpRelEntrConeQuad <- setClass("OpRelEntrConeQuad", representation(X = "ConstValORExpr", Y = "ConstValORExpr", Z = "ConstValORExpr", m = "numeric", k = "numeric"),
                               validity = function(object) {
                                 if(as.integer(object@m) != object@m)
                                   stop("[OpRelEntrConeQuad: m] The argument m must be an integer")
                                 if(object@m <= 0)
                                   stop("[OpRelEntrConeQuad: m] The argument m must be positive")
                                 if(as.integer(object@k) != object@k)
                                   stop("[OpRelEntrConeQuad: k] The argument k must be an integer")
                                 if(object@k <= 0)
                                   stop("[OpRelEntrConeQuad: k] The argument k must be positive")
                                 return(TRUE)
                               }, contains = "Constraint")

#' @param X The variable \eqn{X} in the (approximate) operator relative entropy cone.
#' @param Y The variable \eqn{y} in the (approximate) operator relative entropy cone.
#' @param Z The variable \eqn{z} in the (approximate) operator relative entropy cone.
#' @param m A positive integer that controls the number of quadrature nodes used in a local approximation of the matrix logarithm. Increasing this value results in better local approximations, but does not significantly expand the region of inputs for which the approximation is effective.
#' @param k A positive integer that sets the number of scaling points about which the quadrature approximation is performed. Increasing this value will expand the region of inputs over which the approximation is effective.
#' @param constr_id (Optional) A numeric value representing the constraint ID.
#' @rdname OpRelEntrConeQuad-class
## #' @export
OpRelEntrConeQuad <- function(X, Y, Z, m, k, constr_id = NA_integer_) { .OpRelEntrConeQuad(X = X, Y = Y, Z = Z, m = m, k = k, constr_id = constr_id) }

setMethod("initialize", "OpRelEntrConeQuad", function(.Object, ..., X, Y, Z, m, k) {
  .Object@X <- as.Constant(X)
  .Object@Y <- as.Constant(Y)
  .Object@Z <- as.Constant(Z)
  if(!is_hermitian(X) || !is_hermitian(Y) || !is_hermitian(Z))
    stop(paste(c("One of the input matrices has not explicitly been declared as symmetric or",
                 "Hermitian. If the inputs are Variable objects, try declaring them with the",
                 "symmetric=True or Hermitian=True properties. If the inputs are general ",
                 "Expression objects that are known to be symmetric or Hermitian, then you",
                 "can wrap them with the symmetric_wrap and hermitian_wrap atoms. Failure to",
                 "do one of these things will cause this function to impose a symmetry or",
                 "conjugate-symmetry constraint internally, in a way that is very",
                 "inefficient."), collapse = "\n"))
  .Object@m <- m
  .Object@k <- k
  Xs <- dim(.Object@X)
  Ys <- dim(.Object@Y)
  Zs <- dim(.Object@Z)
  if(!all(Xs == Ys) || !all(Xs == Zs))
    stop("All arguments must have the same shapes. Provided arguments have shapes ", Xs, ", ", Ys, ", and ", Zs)
  callNextMethod(.Object, ..., args = list(.Object@X, .Object@Y, .Object@Z))
})

#' @describeIn OpRelEntrConeQuad Information needed to reconstruct the object aside from the args.
setMethod("get_data", "OpRelEntrConeQuad", function(object) { list(object@m, object@k, id(object)) })

setMethod("show", "OpRelEntrConeQuad", function(object) {
  print(paste("OpRelEntrConeQuad(", as.character(object@X), ", ", as.character(object@Y), ", ", as.character(object@Z), ", ", object@m, ", ", object@k, ")", sep = ""))
})

#' @rdname OpRelEntrConeQuad-class
setMethod("as.character", "OpRelEntrConeQuad", function(x) {
  paste("OpRelEntrConeQuad(", as.character(x@X), ", ", as.character(x@Y), ", ", as.character(x@Z), ", ", x@m, ", ", x@k, ")", sep = "")
})

setMethod("residual", "OpRelEntrConeQuad", function(object) {
  stop("Unimplemented")
})

#' @describeIn OpRelEntrConeQuad The number of entries in the combined cones.
setMethod("size", "OpRelEntrConeQuad", function(object) { 3*num_cones(object) })

#' @describeIn OpRelEntrConeQuad The number of elementwise cones.
setMethod("num_cones", "OpRelEntrConeQuad", function(object) { size(object@X) })

#' @describeIn OpRelEntrConeQuad The dimensions of the exponential cones.
setMethod("cone_sizes", "OpRelEntrConeQuad", function(object) { rep(3, num_cones(object)) })

#' @describeIn OpRelEntrConeQuad An exponential constraint is DCP if each argument is affine.
setMethod("is_dcp", "OpRelEntrConeQuad", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping
    return(all(sapply(object@args, is_affine)))
  }
  return(all(sapply(object@args, is_affine)))
})

#' @describeIn OpRelEntrConeQuad Is the constraint DGP?
setMethod("is_dgp", "OpRelEntrConeQuad", function(object, dpp = FALSE) { FALSE })

#' @describeIn OpRelEntrConeQuad Is the constraint DQCP?
setMethod("is_dqcp", "OpRelEntrConeQuad", function(object) { is_dcp(object) })

setMethod("dim", "OpRelEntrConeQuad", function(x) { c(3, dim(x@X)) })

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn OpRelEntrConeQuad Replaces the dual values of an exponential cone constraint.
setReplaceMethod("dual_value", "OpRelEntrConeQuad", function(object, value) {
  stop("Unimplemented")
})
