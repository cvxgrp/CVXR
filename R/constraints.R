#'
#' The Constraint class.
#'
#' This virtual class represents a mathematical constraint.
#'
#' @name Constraint-class
#' @aliases Constraint
#' @rdname Constraint-class
setClass("Constraint", representation(dual_variables = "list"), prototype(dual_variables = list()), contains = "Canonical")

setMethod("initialize", "Constraint", function(.Object, ..., dual_variables = list()) {
  .Object <- callNextMethod(.Object, ...)
  # .Object@dual_variables <- lapply(.Object@args, function(arg) { Variable(dim(arg)) })
  .Object@dual_variables <- lapply(.Object@args, function(arg) { new("Variable", dim = dim(arg)) })
  return(.Object)
})

#' @param x,object A \linkS4class{Constraint} object.
#' @rdname Constraint-class
setMethod("as.character", "Constraint", function(x) { name(x) })
setMethod("show", "Constraint", function(object) {
  print(paste(class(object), "(", as.character(object@args[[1]]), ")", sep = ""))
})

#' @describeIn Constraint The dimensions of the constrained expression.
setMethod("dim", "Constraint", function(x) { dim(x@args[[1]]) })
#' @describeIn Constraint The size of the constrained expression.
setMethod("size", "Constraint", function(object) { size(object@args[[1]]) })
#' @describeIn Constraint Is the constraint real?
setMethod("is_real", "Constraint", function(object) { !is_complex(object) })
#' @describeIn Constraint Is the constraint imaginary?
setMethod("is_imag", "Constraint", function(object) { all(sapply(object@args, is_imag)) })
#' @describeIn Constraint Is the constraint complex?
setMethod("is_complex", "Constraint", function(object) { any(sapply(object@args, is_complex)) })
#' @describeIn Constraint Is the constraint DCP?
setMethod("is_dcp", "Constraint", function(object) { stop("Unimplemented") })
#' @describeIn Constraint Is the constraint DGP?
setMethod("is_dgp", "Constraint", function(object) { stop("Unimplemented") })
#' @describeIn Constraint The residual of a constraint
setMethod("residual", "Constraint", function(object) { stop("Unimplemented") })

#' @describeIn Constraint The violation of a constraint.
setMethod("violation", "Constraint", function(object) {
  resid <- residual(object)
  if(any(is.na(resid)))
    stop("Cannot compute the violation of a constraint whose expression is NA-valued.")
  return(resid)
})

#' @param tolerance The tolerance for checking if the constraint is violated.
#' @describeIn Constraint The value of a constraint.
setMethod("constr_value", "Constraint", function(object, tolerance = 1e-8) {
  resid <- residual(object)
  if(any(is.na(resid)))
    stop("Cannot compute the value of a constraint whose expression is NA-valued.")
  return(all(resid <= tolerance))
})


#' 
#' A Class Union of List and Constraint
#' 
#' @name ListORConstr-class
#' @rdname ListORConstr-class
setClassUnion("ListORConstr", c("list", "Constraint"))

# Helper function since syntax is different for LinOp (list) vs. Constraint object
#' @param object A list or \linkS4class{Constraint} object.
#' @describeIn ListORConstr Returns the ID associated with the list or constraint.
setMethod("id", "ListORConstr", function(object) {
  if(is.list(object))
    object$constr_id
  else
    object@id
})

#' @describeIn Constraint Information needed to reconstruct the object aside from the args.
setMethod("get_data", "Constraint", function(object) { list(id(object)) })
#' @describeIn Constraint The dual values of a constraint.
setMethod("dual_value", "Constraint", function(object) { value(object@dual_variables[[1]]) })
#' @param value A numeric scalar, vector, or matrix.
#' @describeIn Constraint Replaces the dual values of a constraint..
setReplaceMethod("dual_value", "Constraint", function(object, value) {
  object@dual_variables[[1]] <- value
  object
})

#' 
#' The ZeroConstraint class
#' 
#' @rdname ZeroConstraint-class
.ZeroConstraint <- setClass("ZeroConstraint", representation(expr = "Expression"), contains = "Constraint")
ZeroConstraint <- function(expr, id = NA_integer_) { .ZeroConstraint(expr = expr, id = id) }

setMethod("initialize", "ZeroConstraint", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., args = list(expr))
})

#' @param x,object A \linkS4class{ZeroConstraint} object.
#' @describeIn ZeroConstraint The string representation of the constraint.
setMethod("name", "ZeroConstraint", function(x) {
  # paste(as.character(x@args[[1]]), "== 0")
  paste(name(x@args[[1]]), "== 0")
})

#' @describeIn ZeroConstraint The dimensions of the constrained expression.
setMethod("dim", "ZeroConstraint", function(x) { dim(x@args[[1]]) })
#' @describeIn Constraint The size of the constrained expression.
setMethod("size", "ZeroConstraint", function(object) { size(object@args[[1]]) })
#' @describeIn ZeroConstraint Is the constraint DCP?
setMethod("is_dcp", "ZeroConstraint", function(object) { is_affine(object@args[[1]]) })
#' @describeIn ZeroConstraint Is the constraint DGP?
setMethod("is_dgp", "ZeroConstraint", function(object) { FALSE })
#' @describeIn ZeroConstraint The residual of a constraint
setMethod("residual", "ZeroConstraint", function(object) {
  val <- value(expr(object))
  if(any(is.na(val)))
    return(NA_real_)
  return(abs(val))
})

#' @describeIn ZeroConstraint The graph implementation of the object. 
setMethod("canonicalize", "ZeroConstraint", function(object) {
  canon <- canonical_form(object@args[[1]])
  obj <- canon[[1]]
  constraints <- canon[[2]]
  dual_holder <- create_eq(obj, constr_id = id(object))
  return(list(NA, c(constraints, list(dual_holder))))
})

#' 
#' The EqConstraint class
#' 
#' @rdname EqConstraint-class
.EqConstraint <- setClass("EqConstraint", representation(lhs = "ConstValORExpr", rhs = "ConstValORExpr", expr = "ConstValORExpr"), prototype(expr = NA_real_), contains = "Constraint")
EqConstraint <- function(lhs, rhs, id = NA_integer_) { .EqConstraint(lhs = lhs, rhs = rhs, id = id) }

setMethod("initialize", "EqConstraint", function(.Object, ..., lhs, rhs, expr = NA_real_) {
  .Object@lhs <- lhs
  .Object@rhs <- rhs
  .Object@expr <- lhs - rhs
  callNextMethod(.Object, ..., args = list(lhs, rhs))
})

setMethod(".construct_dual_variables", "EqConstraint", function(object, args) {
  callNextMethod(object, list(object@expr))
})

#' @param x,object A \linkS4class{EqConstraint} object.
#' @describeIn EqConstraint The string representation of the constraint.
setMethod("name", "EqConstraint", function(x) {
  # paste(as.character(x@args[[1]]), "==", as.character(x@args[[2]]))
  paste(name(x@args[[1]]), "==", name(x@args[[2]]))
})

#' @describeIn EqConstraint The dimensions of the constrained expression.
setMethod("dim", "EqConstraint", function(x) { dim(x@expr) })
#' @describeIn EqConstraint The size of the constrained expression.
setMethod("size", "EqConstraint", function(object) { size(object@expr) })
#' @describeIn EqConstraint The expression to constrain.
setMethod("expr", "EqConstraint", function(object) { object@expr })
#' @describeIn EqConstraint Is the constraint DCP?
setMethod("is_dcp", "EqConstraint", function(object) { is_affine(object@expr) })
#' @describeIn EqConstraint Is the constraint DGP?
setMethod("is_dgp", "EqConstraint", function(object) {
  is_log_log_affine(object@args[[1]]) && is_log_log_affine(object@args[[2]])
})

#' @describeIn EqConstraint The residual of the constraint..
setMethod("residual", "EqConstraint", function(object) {
  val <- value(object@expr)
  if(any(is.na(val)))
    return(NA_real_)
  return(abs(val))
})
#' 
#' The NonPosConstraint class
#' 
#' @rdname NonPosConstraint-class
.NonPosConstraint <- setClass("NonPosConstraint", representation(expr = "Expression"), contains = "Constraint")
NonPosConstraint <- function(expr, id = NA_integer_) { .NonPosConstraint(expr = expr, id = id) }

setMethod("initialize", "NonPosConstraint", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., args = list(expr))
})

#' @param x,object A \linkS4class{NonPosConstraint} object.
#' @describeIn NonPosConstraint The string representation of the constraint.
setMethod("name", "NonPosConstraint", function(x) {
  # paste(as.character(x@args[[1]]), "<= 0")
  paste(name(x@args[[1]]), "<= 0")
})

#' @describeIn NonPosConstraint Is the constraint DCP?
setMethod("is_dcp", "NonPosConstraint", function(object) { is_convex(object@args[[1]]) })
#' @describeIn NonPosConstraint Is the constraint DGP?
setMethod("is_dgp", "NonPosConstraint", function(object) { FALSE })
#' @describeIn NonPosConstraint The graph implementation of the object.
setMethod("canonicalize", "NonPosConstraint", function(object) {
  canon <- canonical_form(object@args[[1]])
  obj <- canon[[1]]
  constraints <- canon[[2]]
  dual_holder <- create_leq(obj, constr_id = id(object))
  return(list(NA, c(constraints, list(dual_holder))))
})

#' @describeIn NonPosConstraint The residual of the constraint.
setMethod("residual", "NonPosConstraint", function(object) {
  val <- value(expr(object))
  if(any(is.na(val)))
    return(NA_real_)
  return(pmax(val, 0))
})

#'
#' The IneqConstraint class
#'
#' @rdname IneqConstraint-class
.IneqConstraint <- setClass("IneqConstraint", representation(lhs = "ConstValORExpr", rhs = "ConstValORExpr", expr = "ConstValORExpr"), prototype(expr = NA_real_), contains = "Constraint")

IneqConstraint <- function(lhs, rhs, id = NA_integer_) { .IneqConstraint(lhs = lhs, rhs = rhs, id = id) }

setMethod("initialize", "IneqConstraint", function(.Object, ..., lhs, rhs, expr = NA_real_) {
  .Object@lhs <- lhs
  .Object@rhs <- rhs
  .Object@expr <- lhs - rhs
  if(is_complex(.Object@expr))
    stop("Inequality constraints cannot be complex.")
  callNextMethod(.Object, ..., args = list(lhs, rhs))
})

#' @param x,object A \linkS4class{IneqConstraint} object.
#' @describeIn IneqConstraint The string representation of the constraint.
setMethod("name", "IneqConstraint", function(x) {
  # paste(as.character(x@args[[1]]), "<=", as.character(x@args[[2]]))
  paste(name(x@args[[1]]), "<=", name(x@args[[2]]))
})

#' @describeIn IneqConstraint The dimensions of the constrained expression.
setMethod("dim", "IneqConstraint", function(x) { dim(x@expr) })

#' @describeIn IneqConstraint The size of the constrained expression.
setMethod("size", "IneqConstraint", function(object) { size(object@expr) })

#' @describeIn IneqConstraint The expression to constrain.
setMethod("expr", "IneqConstraint", function(object) { object@expr })

#' @describeIn IneqConstraint A non-positive constraint is DCP if its argument is convex.
setMethod("is_dcp", "IneqConstraint", function(object) { is_convex(object@expr) })

#' @describeIn IneqConstraint Is the constraint DGP?
setMethod("is_dgp", "IneqConstraint", function(object) {
  is_log_log_convex(object@args[[1]]) && is_log_log_concave(object@args[[2]])
})

#' @describeIn IneqConstraint The residual of the constraint.
setMethod("residual", "IneqConstraint", function(object) {
  val <- value(object@expr)
  if(any(is.na(val)))
    return(NA_real_)
  return(pmax(val, 0))
})

# TODO: Do I need the NonlinearConstraint class?
#'
#' The NonlinearConstraint class.
#'
#' This class represents a nonlinear inequality constraint, \eqn{f(x) \leq 0} where \eqn{f} is twice-differentiable.
#'
#' @slot f A nonlinear function.
#' @slot vars_ A list of variables involved in the function.
#' @slot .x_dim (Internal) The dimensions of a column vector with number of elements equal to the total elements in all the variables.
#' @name NonlinearConstraint-class
#' @aliases NonlinearConstraint
#' @rdname NonlinearConstraint-class
.NonlinearConstraint <- setClass("NonlinearConstraint", representation(f = "function", vars_ = "list", .x_dim = "numeric"),
                                 prototype(.x_dim = NULL), contains = "Constraint")

#' @param f A nonlinear function.
#' @param vars_ A list of variables involved in the function.
#' @param id (Optional) An integer representing the unique ID of the contraint.
#' @rdname NonlinearConstraint-class
NonlinearConstraint <- function(f, vars_, id = NA_integer_) { .NonlinearConstraint(f = f, vars_ = vars_, id = id) }

setMethod("initialize", "NonlinearConstraint", function(.Object, ..., f, vars_) {
  .Object@f <- f
  .Object@vars_ <- vars_

  # The dimensions of vars_ in f(vars_)
  sizes <- sapply(.Object@vars_, function(v) { as.integer(prod(dim(v))) })
  .Object@.x_dim <- c(sum(sizes), 1)
  callNextMethod(.Object, ..., args = .Object@vars_)
})

# Add the block to a slice of the matrix.
setMethod("block_add", "NonlinearConstraint", function(object, mat, block, vert_offset, horiz_offset, rows, cols, vert_step = 1, horiz_step = 1) {
  if(is(mat, "sparseMatrix") && is.matrix(block))
    block <- Matrix(block, sparse = TRUE)
  row_seq <- seq(vert_offset, rows + vert_offset, vert_step)
  col_seq <- seq(horiz_offset, cols + horiz_offset, horiz_step)
  mat[row_seq, col_seq] <- mat[row_seq, col_seq] + block
  mat
})

# Place \code{x_0 = f()} in the vector of all variables.
setMethod("place_x0", "NonlinearConstraint", function(object, big_x, var_offsets) {
  tmp <- object@f()
  m <- tmp[[1]]
  x0 <- tmp[[2]]
  offset <- 0
  for(var in object@args) {
    var_dim <- as.integer(prod(dim(var)))
    var_x0 <- x0[offset:(offset + var_dim)]
    big_x <- block_add(object, big_x, var_x0, var_offsets[get_data(var)], 0, var_dim, 1)
    offset <- offset + var_dim
  }
  big_x
})

# Place \code{Df} in the gradient of all functions.
setMethod("place_Df", "NonlinearConstraint", function(object, big_Df, Df, var_offsets, vert_offset) {
  horiz_offset <- 0
  for(var in object@args) {
    var_dim <- as.integer(prod(dim(var)))
    var_Df <- Df[, horiz_offset:(horiz_offset + var_dim)]
    big_Df <- block_add(object, big_Df, var_Df, vert_offset, var_offsets[get_data(var)], num_cones(object), var_dim)
    horiz_offset <- horiz_offset + var_dim
  }
  big_Df
})

# Place \code{H} in the Hessian of all functions.
setMethod("place_H", "NonlinearConstraint", function(object, big_H, H, var_offsets) {
  offset <- 0
  for(var in object@args) {
    var_dim <- as.integer(prod(dim(var)))
    var_H <- H[offset:(offset + var_dim), offset:(offset + var_dim)]
    big_H <- block_add(object, big_H, var_H, var_offsets[get_data(var)], var_offsets[get_data(var)], var_dim, var_dim)
    offset <- offset + var_dim
  }
  big_H
})

# Extract the function variables from the vector \code{x} of all variables.
setMethod("extract_variables", "NonlinearConstraint", function(object, x, var_offsets) {
  local_x <- matrix(0, nrow = object@.x_dim[1], ncol = object@.x_dim[2])
  offset <- 0
  for(var in object@args) {
    var_dim <- as.integer(prod(dim(var)))
    value <- x[var_offsets[get_data(var)]:(var_offsets[get_data(var)] + var_dim)]
    local_x <- block_add(object, local_x, value, offset, 0, var_dim, 1)
    offset <- offset + var_dim
  }
  local_x
})

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
#' @param id (Optional) A numeric value representing the constraint ID.
#' @rdname ExpCone-class
## #' @export
ExpCone <- function(x, y, z, id = NA_integer_) { .ExpCone(x = x, y = y, z = z, id = id) }

setMethod("initialize", "ExpCone", function(.Object, ..., x, y, z) {
  .Object@x <- x
  .Object@y <- y
  .Object@z <- z
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
setMethod("size", "ExpCone", function(object) {
  # TODO: Use size of dual variable(s) instead.
  sum(cone_sizes(object))
})

#' @describeIn ExpCone The number of elementwise cones.
setMethod("num_cones", "ExpCone", function(object) {
  as.integer(prod(dim(object@args[[1]])))
})

#' @describeIn ExpCone The dimensions of the exponential cones.
setMethod("cone_sizes", "ExpCone", function(object) {
  rep(num_cones(object), 3)
})

#' @describeIn ExpCone An exponential constraint is DCP if each argument is affine.
setMethod("is_dcp", "ExpCone", function(object) {
  all(sapply(object@args, is_affine))
})

#' @describeIn ExpCone Is the constraint DGP?
setMethod("is_dgp", "ExpCone", function(object) { FALSE })

#' @describeIn ExpCone Canonicalizes by converting expressions to LinOps.
setMethod("canonicalize", "ExpCone", function(object) {
  arg_objs <- list()
  arg_constr <- list()
  for(arg in object@args) {
    canon <- canonical_form(arg)
    arg_objs <- c(arg_objs, canon[[1]])
    arg_constr <- c(arg_constr, canon[[2]])
  }
  exp_constr <- do.call(ExpCone, arg_objs)
  list(0, c(list(exp_constr), arg_constr))
})

#'
#' The PSDConstraint class.
#'
#' This class represents the positive semidefinite constraint, \eqn{\frac{1}{2}(X + X^T) \succeq 0}, i.e. \eqn{z^T(X + X^T)z \geq 0} for all \eqn{z}.
#'
#' @slot expr An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{X}.
#' @name PSDConstraint-class
#' @aliases PSDConstraint
#' @rdname PSDConstraint-class
.PSDConstraint <- setClass("PSDConstraint", representation(expr = "ConstValORExpr"),
                           validity = function(object) {
                             expr_dim <- dim(object@expr)
                             if(length(expr_dim) != 2 || expr_dim[1] != expr_dim[2])
                               stop("Non-square matrix in positive definite constraint.")
                             return(TRUE)
                           }, contains = "Constraint")

#' @param expr An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{X}.
#' @param id (Optional) A numeric value representing the constraint ID.
#' @rdname PSDConstraint-class
PSDConstraint <- function(expr, id = NA_integer_) { .PSDConstraint(expr = expr, id = id) }

setMethod("initialize", "PSDConstraint", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., args = list(expr))
})

#' @param x,object A \linkS4class{PSDConstraint} object.
#' @describeIn PSDConstraint The string representation of the constraint.
setMethod("name", "PSDConstraint", function(x) {
  # paste(as.character(x@args[[1]]), ">> 0")
  paste(name(x@args[[1]]), ">> 0")
})

#' @describeIn PSDConstraint The constraint is DCP if the left-hand and right-hand expressions are affine.
setMethod("is_dcp", "PSDConstraint", function(object) { is_affine(object@args[[1]]) })

#' @describeIn PSDConstraint Is the constraint DGP?
setMethod("is_dgp", "PSDConstraint", function(object) { FALSE })

#' @describeIn PSDConstraint A \linkS4class{Expression} representing the residual of the constraint.
setMethod("residual", "PSDConstraint", function(object) {
  val <- value(expr(object))
  if(any(is.na(val)))
    return(NA_real_)
  min_eig <- LambdaMin(object@args[[1]] + t(object@args[[1]]))/2
  value(Neg(min_eig))
})

#' @describeIn PSDConstraint The graph implementation of the object. Marks the top level constraint as the \code{dual_holder} so the dual value will be saved to the \linkS4class{PSDConstraint}.
setMethod("canonicalize", "PSDConstraint", function(object) {
  canon <- canonical_form(object@args[[1]])
  obj <- canon[[1]]
  constraints <- canon[[2]]
  dual_holder <- PSDConstraint(obj, id = id(object))
  return(list(NA, c(constraints, list(dual_holder))))
})

setMethod("format_constr", "PSDConstraint", function(object, eq_constr, leq_constr, dims, solver) {
  .format <- function(object) {
    leq_constr <- create_geq(expr(object), constr_id = object@id)
    return(list(leq_constr))
  }
  new_leq_constr <- .format(object)

  # 0 <= A.
  leq_constr <- c(leq_constr, new_leq_constr)
  # Update dims.
  dims[[PSD_DIM]] <- c(dims[[PSD_DIM]], nrow(object))
  list(eq_constr = eq_constr, leq_constr = leq_constr, dims = dims)
})

#'
#' The SOC class.
#'
#' This class represents a second-order cone constraint, i.e. \eqn{\|x\|_2 \leq t}.
#'
#' @slot t The scalar part of the second-order constraint.
#' @slot X A matrix whose rows/columns are each a cone.
#' @slot axis The dimension along which to slice: \code{1} indicates rows, and \code{2} indicates columns. The default is \code{2}.
#' @name SOC-class
#' @aliases SOC
#' @rdname SOC-class
.SOC <- setClass("SOC", representation(t = "ConstValORExpr", X = "ConstValORExpr", axis = "numeric"),
                        prototype(t = NA_real_, X = NA_real_, axis = 2), contains = "Constraint")

#' @param t The scalar part of the second-order constraint.
#' @param X A matrix whose rows/columns are each a cone.
#' @param axis The dimension along which to slice: \code{1} indicates rows, and \code{2} indicates columns. The default is \code{2}.
#' @param id (Optional) A numeric value representing the constraint ID.
#' @rdname SOC-class
## #' @export
SOC <- function(t, X, axis = 2, id = NA_integer_) { .SOC(t = t, X = X, axis = axis, id = id) }

setMethod("initialize", "SOC", function(.Object, ..., t, X, axis = 2) {
  # TODO: Allow imaginary X.
  # if(!(is.null(dim(t)) || length(dim(t)) == 1))
  t_dim <- dim(t)
  if(!(is.null(t_dim) || length(t_dim) == 1 || (length(t_dim) == 2 && t_dim[2] == 1)))
    stop("t must be a scalar or 1-dimensional vector.")
  .Object@t <- t
  .Object@X <- X
  .Object@axis <- axis
  callNextMethod(.Object, ..., args = list(t, X))
})

#' @param x,object A \linkS4class{SOC} object.
#' @rdname SOC-class
setMethod("as.character", "SOC", function(x) {
  paste("SOC(", as.character(x@t), ", ", as.character(x@X), ")", sep = "")
})

#' @describeIn SOC The residual of the second-order constraint.
setMethod("residual", "SOC", function(object) {
  t <- value(object@args[[1]])
  X <- value(object@args[[2]])
  if(is.na(t) || is.na(X))
    return(NA)
  if(object@axis == 2)
    X <- t(X)

  norms <- apply(X, 1, function(row) { norm(row, "2") })
  zero_indices <- which(X <= -t)[1]
  averaged_indices <- which(X >= abs(t))[1]
  X_proj <- as.matrix(X)
  t_proj <- as.matrix(t)
  X_proj[zero_indices] <- 0
  t_proj[zero_indices] <- 0
  avg_coeff <- 0.5*(1 + t/norms)
  X_proj[averaged_indices] <- avg_coeff * X[averaged_indices]
  t_proj[averaged_indices] <- avg_coeff * t[averaged_indices]

  Xt_diff <- cbind(X, t) - cbind(X_proj, t_proj)
  apply(Xt_diff, 1, function(col) { norm(col, "2") })
})

#' @describeIn SOC Information needed to reconstruct the object aside from the args.
setMethod("get_data", "SOC", function(object) { list(object@axis) })

#' @param eq_constr A list of the equality constraints in the canonical problem.
#' @param leq_constr A list of the inequality constraints in the canonical problem.
#' @param dims A list with the dimensions of the conic constraints.
#' @param solver A string representing the solver to be called.
#' @describeIn SOC Format SOC constraints as inequalities for the solver.
setMethod("format_constr", "SOC", function(object, eq_constr, leq_constr, dims, solver) {
  .format <- function(object) {
    list(list(), format_axis(object@args[[1]], object@args[[2]], object@axis))
  }

  leq_constr <- c(leq_constr, .format(object)[[2]])
  dims[[SOC_DIM]] <- c(dims[[SOC_DIM]], cone_sizes(object))
  list(eq_constr = eq_constr, leq_constr = leq_constr, dims = dims)
})

#' @describeIn SOC The number of elementwise cones.
setMethod("num_cones", "SOC", function(object) { prod(dim(object@args[[1]])) })

#' @describeIn SOC The number of entries in the combined cones.
setMethod("size", "SOC", function(object) {
  # TODO: Use size of dual variable(s) instead.
  sum(cone_sizes(object))
})

#' @describeIn SOC The dimensions of the second-order cones.
setMethod("cone_sizes", "SOC", function(object) {
  if(object@axis == 2)   # Collapse columns.
    idx <- 1
  else if(object@axis == 1)   # Collapse rows.
    idx <- 2
  else
    stop("Unimplemented")
  cone_size <- 1 + dim(object@args[[2]])[idx]
  sapply(1:num_cones(object), function(i) { cone_size })
})

#' @describeIn SOC An SOC constraint is DCP if each of its arguments is affine.
setMethod("is_dcp", "SOC", function(object) {
  all(sapply(object@args, function(arg) { is_affine(arg) }))
})

#' @describeIn SOC Is the constraint DGP?
setMethod("is_dgp", "SOC", function(object) { FALSE })

# TODO: Hack.
#' @describeIn SOC The canonicalization of the constraint.
setMethod("canonicalize", "SOC", function(object) {
  canon_t <- canonical_form(object@args[[1]])
  t <- canon_t[[1]]
  t_cons <- canon_t[[2]]

  canon_X <- canonical_form(object@args[[2]])
  X <- canon_X[[1]]
  X_cons <- canon_X[[2]]

  new_soc <- SOC(t, X, object@axis)
  return(list(NA, c(list(new_soc), t_cons, X_cons)))
})

#'
#' The SOCAxis class.
#'
#' This class represents a second-order cone constraint for each row/column.
#' It Assumes \eqn{t} is a vector the same length as \eqn{X}'s rows (columns) for axis == 1 (2).
#'
#' @slot t The scalar part of the second-order constraint.
#' @slot x_elems A list containing \code{X}, a matrix whose rows/columns are each a cone.
#' @slot axis The dimension across which to take the slice: \code{1} indicates rows, and \code{2} indicates columns.
#' @name SOCAxis-class
#' @aliases SOCAxis
#' @rdname SOCAxis-class
.SOCAxis <- setClass("SOCAxis", representation(x_elems = "list"), prototype(x_elems = list()), contains = "SOC")


#' @param t The scalar part of the second-order constraint.
#' @param X A matrix whose rows/columns are each a cone.
#' @param axis The dimension across which to take the slice: \code{1} indicates rows, and \code{2} indicates columns.
#' @param id (Optional) A numeric value representing the constraint ID.
#' @rdname SOCAxis-class
## #' @export
SOCAxis <- function(t, X, axis, id = NA_integer_) { .SOCAxis(t = t, X = X, axis = axis, id = id) }

setMethod("initialize", "SOCAxis", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)
  .Object@x_elems <- list(.Object@X)
  .Object
})

#' @param x,object A \linkS4class{SOCAxis} object.
#' @rdname SOCAxis-class
setMethod("as.character", "SOCAxis", function(x) {
  paste("SOCAxis(", as.character(x@t), ", ", as.character(x@X), ", <", paste(x@axis, collapse = ", "), ">)", sep = "")
})

#' @param eq_constr A list of the equality constraints in the canonical problem.
#' @param leq_constr A list of the inequality constraints in the canonical problem.
#' @param dims A list with the dimensions of the conic constraints.
#' @param solver A string representing the solver to be called.
#' @describeIn SOCAxis Format SOC constraints as inequalities for the solver.
setMethod("format_constr", "SOCAxis", function(object, eq_constr, leq_constr, dims, solver) {
 .format <- function(object) {
   list(list(), format_axis(object@t, object@x_elems[[1]], object@axis))
 }

 leq_constr <- c(leq_constr, .format(object)[[2]])

 # Update dims
 for(cone_size in size(object))
   dims[[SOC_DIM]] <- c(dims[[SOC_DIM]], cone_size[1])
 list(eq_constr = eq_constr, leq_constr = leq_constr, dims = dims)
})

#' @describeIn SOCAxis The number of elementwise cones.
setMethod("num_cones", "SOCAxis", function(object) { nrow(object@t) })

#' @describeIn SOCAxis The dimensions of a single cone.
setMethod("cone_sizes", "SOCAxis", function(object) {
  if(object@axis == 1)   # Return ncols if applying along each row
    c(1 + ncol(object@x_elems[[1]]), 1)
  else if(object@axis == 2)   # Return nrows if applying along each column
    c(1 + nrow(object@x_elems[[1]]), 1)
  else
    stop("Invalid axis ", object@axis)
})

#' @describeIn SOCAxis The dimensions of the (elementwise) second-order cones.
setMethod("size", "SOCAxis", function(object) {
  cone_size <- cone_size(object)
  lapply(1:num_cones(object), function(i) { cone_size })
})
