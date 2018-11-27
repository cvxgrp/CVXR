#'
#' The Constraint class.
#'
#' This virtual class represents a mathematical constraint.
#'
#' @slot constr_id (Internal) A unique integer identification number used internally.
#' @name Constraint-class
#' @aliases Constraint
#' @rdname Constraint-class
setClass("Constraint", representation(args = "list", constr_id = "integer", dual_variables = "list"),
                       prototype(constr_id = NA_integer_, dual_variables = list()), contains = "Canonical")

setMethod("initialize", "Constraint", function(.Object, ..., args, constr_id = get_id(), dual_variables = list()) {
  .Object@args <- args
  .Object@constr_id <- constr_id
  .Object@dual_variables <- lapply(args, function(arg) { Variable(shape(arg)) })
  callNextMethod(.Object, ...)
})

setMethod("as.character", "Constraint", function(x) { name(x) })
setMethod("show", "Constraint", function(object) {
  paste(class(object), "(", as.character(object@args[[1]]), ")")
})

setMethod("shape", "Constraint", function(object) { shape(object@args[[1]]) })
setMethod("size", "Constraint", function(object) { size(object@args[[1]]) })
setMethod("is_real", "Constraint", function(object) { !is_complex(object) })
setMethod("is_imag", "Constraint", function(object) { all(sapply(object@args, is_imag)) })
setMethod("is_complex", "Constraint", function(object) { any(sapply(object@args, is_complex)) })
setMethod("is_dcp", "Constraint", function(object) { stop("Unimplemented") })
setMethod("residual", "Constraint", function(object) { stop("Unimplemented") })

setMethod("violation", "Constraint", function(object) {
  residual <- object@residual
  if(is.na(residual))
    stop("Cannot compute the violation of a constraint whose expression is NA-valued.")
  return(residual)
})

setMethod("value", "Constraint", function(object, tolerance = 1e-8) {
  residual <- object@residual
  if(is.na(residual))
    stop("Cannot compute the value of a constraint whose expression is NA-valued.")
  return(all(residual <= tolerance))
})

setClassUnion("ListORConstr", c("list", "Constraint"))

# Helper function since syntax is different for LinOp (list) vs. Constraint object
setMethod("id", "ListORConstr", function(object) {
  if(is.list(object))
    object$constr_id
  else
    object@constr_id
})

setMethod("get_data", "Constraint", function(object) { list(id(object)) })
setMethod("dual_value", "Constraint", function(object) { value(object@dual_variables[[1]]) })
setMethod("save_dual_value", "Constraint", function(object, value) {
  object@dual_variables[[1]] <- value
  object
})

.ZeroConstraint <- setClass("ZeroConstraint", representation(expr = "ConstValORExpr"), contains = "Constraint")
ZeroConstraint <- function(expr, constr_id = NA_integer_) { .ZeroConstraint(expr = expr, constr_id = constr_id) }

setMethod("initialize", "ZeroConstraint", function(.Object, ..., expr) {
  callNextMethod(.Object, ..., args = list(expr))
})

setMethod("name", "ZeroConstraint", function(x) {
  paste(as.character(x@args[[1]]), "== 0")
})

setMethod("is_dcp", "ZeroConstraint", function(object) { is_affine(object@args[[1]]) })
setMethod("residual", "ZeroConstraint", function(object) {
  if(is.na(value(object@expr)))
    return(NA_real_)
  return(abs(value(object@expr)))
})

setMethod("canonicalize", "ZeroConstraint", function(object) {
  canon <- canonical_form(object@args[[1]])
  obj <- canon[[1]]
  constraints <- canon[[2]]
  dual_holder <- create_eq(obj, constr_id = id(object))
  return(list(NA, c(constraints, list(dual_holder))))
})

.NonPosConstraint <- setClass("NonPosConstraint", representation(expr = "ConstValORExpr"),
                              validity = function(object) {
                                if(is_complex(object@expr))
                                  stop("Inequality constraints cannot be complex.")
                                return(TRUE)
                              }, contains = "Constraint")
NonPosConstraint <- function(expr, constr_id = NA_integer_) { .NonPosConstraint(expr = expr, constr_id = constr_id) }

setMethod("initialize", "NonPosConstraint", function(.Object, ..., expr) {
  callNextMethod(.Object, ..., args = list(expr))
})

setMethod("name", "NonPosConstraint", function(x) {
  paste(as.character(x@args[[1]]), "<= 0")
})

setMethod("is_dcp", "NonPosConstraint", function(object) { is_convex(object@args[[1]]) })
setMethod("canonicalize", "NonPosConstraint", function(object) {
  canon <- canonical_form(object@args[[1]])
  obj <- canon[[1]]
  constraints <- canon[[2]]
  dual_holder <- create_leq(obj, constr_id = id(object))
  return(list(NA, c(constraints, list(dual_holder))))
})

setMethod("residual", "NonPosConstraint", function(object) {
  if(is.na(value(object@expr)))
    return(NA_real_)
  return(pmax(value(object@expr), 0))
})

# TODO: Do I need the NonlinearConstraint class?
#'
#' The NonlinearConstraint class.
#'
#' This class represents a nonlinear inequality constraint, \eqn{f(x) \leq 0} where \eqn{f} is twice-differentiable.
#'
#' @slot constr_id (Internal) A unique integer identification number used internally.
#' @slot f A nonlinear function.
#' @slot vars_ A list of variables involved in the function.
#' @slot .x_size (Internal) The dimensions of a column vector with number of elements equal to the total elements in all the variables.
#' @name NonlinearConstraint-class
#' @aliases NonlinearConstraint
#' @rdname NonlinearConstraint-class
.NonlinearConstraint <- setClass("NonlinearConstraint", representation(f = "function", vars_ = "list", .x_shape = "numeric"),
                                 prototype(.x_shape = NA_integer_), contains = "Constraint")

#' @param f A nonlinear function.
#' @param vars_ A list of variables involved in the function.
#' @rdname NonlinearConstraint-class
NonlinearConstraint <- function(f, vars_, constr_id = NA_integer_) { .NonlinearConstraint(f = f, vars_ = vars_, constr_id = constr_id) }

setMethod("initialize", "NonlinearConstraint", function(.Object, ..., f, vars_) {
  .Object@f <- f
  .Object@vars_ <- vars_

  # The shape of vars_ in f(vars_)
  sizes <- sapply(.Object@vars_, function(v) { as.integer(prod(shape(v))) })
  .Object@x_shape <- c(sum(sizes), 1)
  callNextMethod(.Object, ..., args = .Object@vars_)
})

# Add the block to a slice of the matrix.
setMethod("block_add", "NonlinearConstraint", function(object, matrix, block, vert_offset, horiz_offset, rows, cols, vert_step = 1, horiz_step = 1) {
  if(is(matrix, "sparseMatrix") && is.matrix(block))
    block <- Matrix(block, sparse = TRUE)
  row_seq <- seq(vert_offset, rows + vert_offset, vert_step)
  col_seq <- seq(horiz_offset, cols + horiz_offset, horiz_step)
  matrix[row_seq, col_seq] <- matrix[row_seq, col_seq] + block
  matrix
})

# Place \code{x_0 = f()} in the vector of all variables.
setMethod("place_x0", "NonlinearConstraint", function(object, big_x, var_offsets) {
  tmp <- object@f()
  m <- tmp[[1]]
  x0 <- tmp[[2]]
  offset <- 0
  for(var in object@args) {
    var_shape <- as.integer(prod(shape(var)))
    var_x0 <- x0[offset:(offset + var_shape)]
    big_x <- block_add(object, big_x, var_x0, var_offsets[get_data(var)], 0, var_shape, 1)
    offset <- offset + var_shape
  }
  big_x
})

# Place \code{Df} in the gradient of all functions.
setMethod("place_Df", "NonlinearConstraint", function(object, big_Df, Df, var_offsets, vert_offset) {
  horiz_offset <- 0
  for(var in object@args) {
    var_shape <- as.integer(prod(shape(var)))
    var_Df <- Df[, horiz_offset:(horiz_offset + var_shape)]
    big_Df <- block_add(object, big_Df, var_Df, vert_offset, var_offsets[get_data(var)], num_cones(object), var_shape)
    horiz_offset <- horiz_offset + var_shape
  }
  big_Df
})

# Place \code{H} in the Hessian of all functions.
setMethod("place_H", "NonlinearConstraint", function(object, big_H, H, var_offsets) {
  offset <- 0
  for(var in object@args) {
    var_size <- as.integer(prod(shape(var)))
    var_H <- H[offset:(offset + var_shape), offset:(offset + var_shape)]
    big_H <- block_add(object, big_H, var_H, var_offsets[get_data(var)], var_offsets[get_data(var)], var_shape, var_shape)
    offset <- offset + var_shape
  }
  big_H
})

# Extract the function variables from the vector \code{x} of all variables.
setMethod("extract_variables", "NonlinearConstraint", function(object, x, var_offsets) {
  local_x <- matrix(0, nrow = object@.x_shape[1], ncol = object@.x_shape[2])
  offset <- 0
  for(var in object@args) {
    var_shape <- as.integer(prod(shape(var)))
    value <- x[var_offsets[get_data(var)]:(var_offsets[get_data(var)] + var_shape)]
    local_x <- block_add(object, local_x, value, offset, 0, var_shape, 1)
    offset <- offset + var_shape
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
#' @slot constr_id (Internal) A unique integer identification number used internally.
#' @slot x The variable \eqn{x} in the exponential cone.
#' @slot y The variable \eqn{y} in the exponential cone.
#' @slot z The variable \eqn{z} in the exponential cone.
#' @name ExpCone-class
#' @aliases ExpCone
#' @rdname ExpCone-class
.ExpCone <- setClass("ExpCone", representation(x = "list", y = "list", z = "list"), contains = "NonlinearConstraint")

#' @param x The variable \eqn{x} in the exponential cone.
#' @param y The variable \eqn{y} in the exponential cone.
#' @param z The variable \eqn{z} in the exponential cone.
#' @rdname ExpCone-class
## #' @export
ExpCone <- function(x, y, z) { .ExpCone(x = x, y = y, z = z) }

setMethod("initialize", "ExpCone", function(.Object, ..., x, y, z) {
  .Object@x <- x
  .Object@y <- y
  .Object@z <- z
  callNextMethod(.Object, ..., f = solver_hook(.Object), vars_ = list(.Object@x, .Object@y, .Object@z))
})

setMethod("show", "ExpCone", function(object) {
  paste("ExpCone(", as.character(x@a), ", ", as.character(x@b), ", ", as.character(x@c), ")", sep = "")
})

#' @rdname ExpCone-class
setMethod("as.character", "ExpCone", function(x) {
  paste("ExpCone(", as.character(x@a), ", ", as.character(x@b), ", ", as.character(x@c), ")", sep = "")
})

#' @param x,object A \linkS4class{ExpCone} object.
#' @describeIn ExpCone The size of the \code{x} argument.
setMethod("residual", "ExpCone", function(object) {
  # TODO: The projection should be implemented directly.
  if(is.na(value(object@x)) || is.na(value(object@y)) || is.na(value(object@z)))
    return(NA_real_)
  x <- Variable(shape(object@x))
  y <- Variable(shape(object@y))
  z <- Variable(shape(object@z))
  constr <- list(ExpCone(x, y, z))
  obj <- Minimize(norm2(hstack(list(x, y, z)) - hstack(list(value(object@x), value(object@y), value(object@z)))))
  prob <- Problem(obj, constr)
  return(solve(prob))
})

#' @param eq_constr A list of the equality constraints in the canonical problem.
#' @param leq_constr A list of the inequality constraints in the canonical problem.
#' @param dims A list with the dimensions of the conic constraints.
#' @param solver A string representing the solver to be called.
#' @describeIn ExpCone Format exponential cone constraints for the solver.
setMethod("format_constr", "ExpCone", function(object, eq_constr, leq_constr, dims, solver) {
  .cvxopt_format <- function(object) {
    constraints <- list()
    for(i in 1:length(object@vars_)) {
      var <- object@vars_[[i]]
      if(var$type != VARIABLE) {
        lone_var <- create_var(var$size)
        constraints <- c(constraints, list(create_eq(lone_var, var)))
        object@vars_[[i]] <- lone_var
      }
    }
    list(constraints, list())
  }

  .ecos_format <- function(object) {
    list(list(), format_elemwise(list(object@x, object@y, object@z)))
  }

  .scs_format <- function(object) {
    list(list(), format_elemwise(list(object@z, object@y, object@z)))
  }

  if(is(solver, "CVXOPT"))
    stop("CVXOPT formatting has not been implemented")
    # tmp <- .cvxopt_format(object)
    # object <- tmp[[3]]   # TODO: Need to pass along updated constraint object
    # eq_constr <- c(eq_constr, tmp[[1]])
  else if(is(solver, "SCS") || is(solver, "JULIA_OPT"))
    leq_constr <- c(leq_constr, .scs_format(object)[[2]])
  else if(is(solver, "ECOS"))
    leq_constr <- c(leq_constr, .ecos_format(object)[[2]])
  else
    stop("Solver does not support exponential cone")
  # Update dims
  dims[[EXP_DIM]] <- dims[[EXP_DIM]] + num_cones(object)
  list(eq_constr = eq_constr, leq_constr = leq_constr, dims = dims)
})

#' @describeIn ExpCone The number of entries in the combined cones.
setMethod("size", "ExpCone", function(object) {
  # TODO: Use size of dual variable(s) instead.
  sum(cone_sizes(object))
})

#' @describeIn ExpCone The number of elementwise cones.
setMethod("num_cones", "ExpCone", function(object) {
  as.integer(prod(shape(object@args[[1]])))
})

#' @describeIn ExpCone The dimensions of the exponential cones.
setMethod("cone_sizes", "ExpCone", function(object) {
  rep(num_cones(object), 3)
})

#' @describeIn ExpCone An exponential constraint is DCP if each argument is affine.
setMethod("is_dcp", "ExpCone", function(object) {
  all(sapply(object@args, function(arg) { is_affine(arg) }))
})

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

# A function used by CVXOPT's nonlinear solver.
# Based on f(x,y,z) = y * log(y) + x - y * log(z)
setMethod("solver_hook", "ExpCone", function(object, vars_ = NA, scaling = NA) {
  entries <- as.integer(prod(shape(object)))
  if(is.na(vars_)) {
    x_init <- rep(0, entries)
    y_init <- rep(0.5, entries)
    z_init <- rep(1.0, entries)
    return(list(entries, matrix(x_init + y_init + z_init)))
  }

  # Unpack vars_
  x <- vars_[1:entries]
  y <- vars_[entries:(2*entries)]
  z <- vars_[(2*entries):length(vars_)]

  # Out of domain.
  # TODO: What if y == 0.0?
  if(min(y) <= 0.0 || min(z) <= 0.0)
    return(NA)

  # Evaluate the function.
  f <- matrix(0, nrow = entries, ncol = 1)
  for(i in 1:entries)
    f[i] <- x[i] - y[i]*log(z[i]) + y[i]*log(y[i])

  # Compute the gradient.
  Df <- matrix(0, nrow = entries, ncol = 3*entries)
  for(i in 1:entries) {
    Df[i,i] <- 1.0
    Df[i, entries+i] <- log(y[i]) - log(z[i]) + 1.0
    Df[i, 2*entries+i] <- -y[i]/z[i]
  }

  if(is.na(scaling))
    return(list(f, Df))

  # Compute the Hessian.
  big_H <- sparseMatrix(i = c(), j = c(), dims = c(3*entries, 3*entries))
  for(i in 1:entries) {
    H <- rbind(c(0.0, 0.0, 0.0),
               c(0.0, 1.0/y[i], -1.0/z[i]),
               c(0.0, -1.0/z[i], y[i]/(z[i]^2)))
    idx <- seq(i, 3*entries, entries)
    big_H[idx, idx] <- scaling[i]*H
  }
  return(list(f, Df, big_H))
})

#'
#' The PSDConstraint class.
#'
#' This class represents the positive semidefinite constraint, \eqn{X \succeq Y}, i.e. \eqn{z^T(X - Y)z \geq 0} for all \eqn{z}.
#'
#' @slot constr_id (Internal) A unique integer identification number used internally.
#' @slot lh_exp An \linkS4class{Expression}, numeric element, vector, or matrix representing the left-hand side of the inequality.
#' @slot rh_exp An \linkS4class{Expression}, numeric element, vector, or matrix representing the right-hand side of the inequality.
#' @slot args (Internal) A list that holds \code{lh_exp} and \code{rh_exp} for internal use.
#' @slot .expr (Internal) An \linkS4class{Expression} representing \code{lh_exp - rh_exp} for internal use.
#' @slot dual_variable (Internal) A \linkS4class{Variable} representing the dual variable associated with the constraint.
#' @name PSDConstraint-class
#' @aliases PSDConstraint
#' @rdname PSDConstraint-class
.PSDConstraint <- setClass("PSDConstraint", representation(expr = "ConstValORExpr"),
                           validity = function(object) {
                             shape <- shape(object@expr)
                             if(length(shape) != 2 || shape[1] != shape[2])
                               stop("Non-square matrix in positive definite constraint.")
                             return(TRUE)
                           }, contains = "Constraint")

#' @param lh_exp An \linkS4class{Expression}, numeric element, vector, or matrix representing the left-hand side of the inequality.
#' @param rh_exp An \linkS4class{Expression}, numeric element, vector, or matrix representing the right-hand side of the inequality.
#' @rdname PSDConstraint-class
PSDConstraint <- function(expr, constr_id = NA_integer_) { .PSDConstraint(expr = expr, constr_id = constr_id) }

setMethod("name", "PSDConstraint", function(x) {
  paste(as.character(x@args[[1]]), ">> 0")
})

#' @param object A \linkS4class{PSDConstraint} object.
#' @describeIn PSDConstraint The constraint is DCP if the left-hand and right-hand expressions are affine.
setMethod("is_dcp", "PSDConstraint", function(object) { is_affine(object@args[[1]]) })

#' @describeIn PSDConstraint A \linkS4class{Expression} representing the residual of the constraint.
setMethod("residual", "PSDConstraint", function(object) {
  if(is.na(value(object@expr)))
    return(NA_real_)
  min_eig <- LambdaMin(object@args[[1]] + t(object@args[[1]]))/2
  value(Neg(min_eig))
})

#' @describeIn PSDConstraint The graph implementation of the object. Marks the top level constraint as the \code{dual_holder} so the dual value will be saved to the \linkS4class{PSDConstraint}.
setMethod("canonicalize", "PSDConstraint", function(object) {
  canon <- canonical_form(object@args[[1]])
  obj <- canon[[1]]
  constraints <- canon[[2]]
  dual_holder <- PSDConstraint(obj, constr_id = id(object))
  return(list(NA, c(constraints, list(dual_holder))))
})

setMethod("format_constr", "PSDConstraint", function(object, eq_constr, leq_constr, dims, solver) {
  .format <- function(object) {
    leq_constr <- create_geq(object@expr, constr_id = object@constr_id)
    return(list(leq_constr))
  }
  new_leq_constr <- .format(object)

  # 0 <= A.
  leq_constr <- c(leq_constr, new_leq_constr)
  # Update dims.
  dims[[PSD_DIM]] <- c(dims[[PSD_DIM]], shape(object)[1])
  list(eq_constr = eq_constr, leq_constr = leq_constr, dims = dims)
})

#'
#' The SOC class.
#'
#' This class represents a second-order cone constraint, i.e. \eqn{\|x\|_2 \leq t}.
#'
#' @slot constr_id (Internal) A unique integer identification number used internally.
#' @slot t The scalar part of the second-order constraint.
#' @slot x_elems A list containing the elements of the vector part of the constraint.
#' @name SOC-class
#' @aliases SOC
#' @rdname SOC-class
.SOC <- setClass("SOC", representation(t = "ConstValORExpr", X = "ConstValORExpr", axis = "numeric"),
                        prototype(t = NA_real_, X = NA_real_, axis = 2), contains = "Constraint")

#' @param t The scalar part of the second-order constraint.
#' @param x_elems A list containing the elements of the vector part of the constraint.
#' @rdname SOC-class
## #' @export
SOC <- function(t, X, axis = 2, constr_id = NA_integer_) { .SOC(t = t, X = X, axis = axis, constr_id = constr_id) }

setMethod("initialize", "SOC", function(.Object, ..., t, X, axis = 2) {
  # TODO: Allow imaginary X.
  if(!(is.null(shape(t)) || length(shape(t)) == 1))
    stop("t must be a scalar or 1-dimensional vector.")
  .Object@t <- t
  .Object@X <- X
  .Object@axis <- axis
  callNextMethod(.Object, ..., args = list(t, X))
})

#' @param x,object A \linkS4class{SOC} object.
#' @rdname SOC-class
setMethod("as.character", "SOC", function(x) {
  paste("SOC(", as.character(x@t), ", <", paste(lapply(x@x_elems, function(elem) { as.character(elem) }), collapse = ", "), ">)", sep = "")
})

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
setMethod("num_cones", "SOC", function(object) { prod(shape(object@args[[1]])) })

#' @describeIn SOC The number of entries in the combined cones.
setMethod("size", "SOC", function(object) {
  # TODO: Use size of dual variable(s) instead.
  sum(cone_sizes(object))
})

#' @describeIn SOC The dimensions of the second-order cones.
setMethod("cone_sizes", "SOC", function(object) {
  if(object@axis == 2)   # Collapse columns.
    idx <- 1
  else   # Collapse rows.
    idx <- 2
  if(object@axis == 2)
    cone_size <- 1 + shape(object@args[[2]])[1]
  else if(object@axis == 1)
    cone_size <- 1 + shape(object@args[[2]])[2]
  else
    stop("Invalid axis ", object@axis)
  lapply(1:num_cones(object), function(i) { cone_size })
})

#' @describeIn SOC An SOC constraint is DCP if each of its arguments is affine.
setMethod("is_dcp", "SOC", function(object) {
  all(sapply(object@args, function(arg) { is_affine(arg) }))
})

# TODO: Hack.
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
#' @slot constr_id (Internal) A unique integer identification number used internally.
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
#' @rdname SOCAxis-class
## #' @export
SOCAxis <- function(t, X, axis) { .SOCAxis(t = t, X = X, axis = axis) }

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
setMethod("num_cones", "SOCAxis", function(object) { shape(object@t)[1] })

#' @describeIn SOCAxis The dimensions of a single cone.
setMethod("cone_size", "SOCAxis", function(object) {
  if(object@axis == 1)   # Return ncols if applying along each row
    c(1 + shape(object@x_elems[[1]])[2], 1)
  else if(object@axis == 2)   # Return nrows if applying along each column
    c(1 + shape(object@x_elems[[1]])[1], 1)
  else
    stop("Invalid axis ", object@axis)
})

#' @describeIn SOCAxis The dimensions of the (elementwise) second-order cones.
setMethod("size", "SOCAxis", function(object) {
  cone_size <- cone_size(object)
  lapply(1:num_cones(object), function(i) { cone_size })
})
