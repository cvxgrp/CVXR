#'
#' The Constraint class.
#'
#' This virtual class represents a mathematical constraint.
#'
#' @name Constraint-class
#' @aliases Constraint
#' @rdname Constraint-class
setClass("Constraint", representation(dual_variables = "list", constr_id = "integer"), prototype(dual_variables = list(), constr_id = NA_integer_), contains = "Canonical")

setMethod("initialize", "Constraint", function(.Object, ..., dual_variables = list(), constr_id = NA_integer_) {
  .Object <- callNextMethod(.Object, ...)
  # TODO: Cast constants
  # .Object@args <- lapply(args, as.Constant)
  .Object@constr_id <- ifelse(is.na(constr_id), get_id(), constr_id)
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

#' @param dpp A logical value indicating whether we are solving a disciplined parameterized program (DPP).
#' @describeIn Constraint Is the constraint DCP?
setMethod("is_dcp", "Constraint", function(object, dpp = FALSE) { stop("Unimplemented") })

#' @describeIn Constraint Is the constraint DGP?
setMethod("is_dgp", "Constraint", function(object, dpp = FALSE) { stop("Unimplemented") })

#' @param context Must be either 'dcp' (disciplined convex program) or 'dgp' (disciplined geometric program).
#' @describeIn Constraint is the constraint DPP in the given context?
setMethod("is_dpp", "Constraint", function(object, context = "dcp") {
  if(tolower(context) == "dcp")
    return(is_dcp(object, dpp = TRUE))
  else if(tolower(context) == "dgp")
    return(is_dpp(object, dpp = TRUE))
  else
    stop("Unsupported context ", context)
})

#' @describeIn Constraint The residual of a constraint
setMethod("residual", "Constraint", function(object) { stop("Unimplemented") })

#' @describeIn Constraint The numeric residual of a constraint. The violation is defined as the 
#' distance between the constrained expression's value and its projection onto the domain of the
#' constraint: ||\\Pi(v) - v||_2^2, where `v` is the value of the constrained expression and
#' `\\Pi` is the projection operator onto the constraint's domain.
setMethod("violation", "Constraint", function(object) {
  resid <- residual(object)
  if(any(is.na(resid)))
    stop("Cannot compute the violation of a constraint whose expression is NA-valued.")
  return(resid)
})

#' @param tolerance The tolerance for checking if the constraint is violated.
#' @describeIn Constraint Checks whether the constraint violation is less than a tolerance.
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
    object@constr_id
})

#' @describeIn Constraint Information needed to reconstruct the object aside from the args.
setMethod("get_data", "Constraint", function(object) { list(id(object)) })

#' @describeIn Constraint The dual values of a constraint.
setMethod("dual_value", "Constraint", function(object) {
  dual_vals <- lapply(object@dual_variables, value)
  if(length(dual_vals) == 1)
    return(dual_vals[[1]])
  else
    return(dual_vals)
})

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn Constraint Replaces the dual values of a constraint.
setReplaceMethod("dual_value", "Constraint", function(object, value) {
  value(object@dual_variables[[1]]) <- value
  return(object)
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
#' A constraint of the form \eqn{x \leq 0}.
#' 
#' The preferred way of creating a NonPosConstraint constraint is through
#' operator overloading. To constrain an expression x to be non-positive,
#' simply write \code{x <= 0}; to constrain x to be non-negative, write
#' \code{x >= 0}. The former creates a NonPosConstraint constraint with x
#' as its argument, while the latter creates one with -x as its argument.
#' Strict inequalities are not supported, as they do not make sense in a
#' numerical setting.
#' 
#' @rdname NonPosConstraint-class
.NonPosConstraint <- setClass("NonPosConstraint", representation(expr = "Expression"), contains = "Constraint")
NonPosConstraint <- function(expr, constr_id = NA_integer_) { .NonPosConstraint(expr = expr, constr_id = constr_id) }

setMethod("initialize", "NonPosConstraint", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., args = list(expr))
  if(!is_real(.Object@args[[1]]))
    stop("Input to NonPosConstraint must be real")
})

#' @param x,object A \linkS4class{NonPosConstraint} object.
#' @describeIn NonPosConstraint The string representation of the constraint.
setMethod("name", "NonPosConstraint", function(x) {
  # paste(as.character(x@args[[1]]), "<= 0")
  paste(name(x@args[[1]]), "<= 0")
})

#' @describeIn NonPosConstraint A non-positive constraint is DCP if its argument is convex.
setMethod("is_dcp", "NonPosConstraint", function(object, dpp = FALSE) { 
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(is_convex(object@args[[1]]))
  }
  return(is_convex(object@args[[1]]))
})

#' @describeIn NonPosConstraint Is the constraint DGP?
setMethod("is_dgp", "NonPosConstraint", function(object, dpp = FALSE) { FALSE })

#' @describeIn NonPosConstraint Is the constraint DQCP?
setMethod("is_dqcp", "NonPosConstraint", function(object) { is_quasiconvex(object@args[[1]]) })

#' @describeIn NonPosConstraint The residual of the constraint.
setMethod("residual", "NonPosConstraint", function(object) {
  val <- value(expr(object))
  if(any(is.na(val)))
    return(NA_real_)
  return(pmax(val, 0))
})

#' @describeIn NonPosConstraint The violation of the constraint.
setMethod("violation", "NonPosConstraint", function(object) {
  resid <- residual(object)
  if(any(is.na(resid)))
    stop("Cannot compute the violation of a constraint whose expression is NA-valued.")
  viol <- base::norm(resid, type = "2")
  return(viol)
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

#'
#' The FiniteSet class.
#' 
#' This class represents a constraint that each entry of an Expression to take a value in a given set of real numbers.
#' 
#' @slot expre The given expression to be constrained. This Expression must be affine.
#' If expre has multiple elements, then the constraint is applied separately to
#' each element, i.e., after solving a problem with this constraint, we should have:
#' \code{for(e in flatten(expre)) { print(value(e) %in% vec) # => TRUE }
#' @slot vec The finite collection of values to which each entry of expre is to be constrained.
#' @slot ineq_form A logical value controlling how this constraint is canonicalized into mixed-integer linear constraints.
#' If TRUE, then we use a formulation with (size(vec) - 1) inequality constraints,
#' one equality constraint, and (size(vec) - 1) binary variables for each element
#' of expre. If FALSE, then we use a formulation with size(vec) binary variables and two
#' equality constraints for each element of expre. Defaults to FALSE. The case ineq_form = TRUE may speed up some mixed-integer
#' solvers that use simple branch and bound methods.
#' @name FiniteSet-class
#' @aliases FiniteSet
#' @rdname FiniteSet-class
.FiniteSet <- setClass("FiniteSet", representation(expre = "ConstValORExpr", vec = "list", ineq_form = "logical"),
                       prototype(ineq_form = FALSE), contains = "Constraint")

#' @param expre An affine Expression object.
#' @param vec The finite collection of values to which each entry of expre is to be constrained.
#' @param ineq_form A logical value controlling how this constraint is canonicalized.
#' @param constr_id (Optional) An integer representing the unique ID of the constraint.
#' @rdname FiniteSet-class
FiniteSet <- function(expre, vec, ineq_form = FALSE, constr_id = NA_integer_) { .FiniteSet(expre = expre, vec = vec, ineq_form = ineq_form, constr_id = constr_id) }

setMethod("initialize", "FiniteSet", function(.Object, ..., expre, vec, ineq_form = FALSE) {
  vec <- flatten(as.Constant(vec))
  if(!is_affine(expre))
    stop("Provided Expression must be affine, but had curvature ", curvature(expre))
    
  # Note: we use the term "expre" rather than "expr" since
  # "expr" is already a property used by all Constraint classes.
  .Object@expre <- expre
  .Object@vec <- vec
  .Object@ineq_form <- ineq_form
  callNextMethod(.Object, ..., args = list(expre, vec))
})

#' @param x,object A \linkS4class{FiniteSet} object.
#' @describeIn FiniteSet The string representation of the constraint.
setMethod("name", "FiniteSet", function(x) {
  paste("FiniteSet(", as.character(x@args[[1]]), ", ", as.character(x@args[[2]]), ")", sep = "")
})

#' @describeIn FiniteSet Information needed to reconstruct the object aside from the args.
setMethod("get_data", "FiniteSet", function(object) { list(object@ineq_form, id(object)) })

#' @describeIn FiniteSet The constraint is DCP if the constrained expression is affine.
setMethod("is_dcp", "FiniteSet", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping
    return(is_affine(object@args[[1]]))
  }
  return(is_affine(object@args[[1]]))
})

#' @describeIn FiniteSet Is the constraint DGP?
setMethod("is_dgp", "FiniteSet", function(object, dpp = FALSE) { FALSE })

#' @describeIn FiniteSet Is the constraint DQCP?
setMethod("is_dqcp", "FiniteSet", function(object) { is_dcp(object) })

#' @describeIn FiniteSet The size of the constrained expression.
setMethod("size", "FiniteSet", function(object) { size(object@expre) })

setMethod("dim", "FiniteSet", function(x) { dim(object@expre) })

#' @describeIn FiniteSet The residual of the constraint.
setMethod("residual", "FiniteSet", function(object) {
  expr_val <- as.vector(value(object@expre))
  vec_val <- value(object@vec)
  resids <- sapply(expr_val, min(abs(val - vec_val)))
  res <- max(resids)
  return(res)
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
#' @param id (Optional) An integer representing the unique ID of the constraint.
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
  # Python code: value = np.reshape(value, newshape=(-1, 3))
  value <- matrix(t(value), ncol = 3, byrow = TRUE)
  
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
    stop("One of the input matrices has not explicitly been declared as symmetric or"
         "Hermitian. If the inputs are Variable objects, try declaring them with the"
         "symmetric=True or Hermitian=True properties. If the inputs are general "
         "Expression objects that are known to be symmetric or Hermitian, then you"
         "can wrap them with the symmetric_wrap and hermitian_wrap atoms. Failure to"
         "do one of these things will cause this function to impose a symmetry or"
         "conjugate-symmetry constraint internally, in a way that is very"
         "inefficient.")
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
