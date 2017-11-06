#'
#' The Constraint class.
#'
#' This virtual class represents a mathematical constraint.
#'
#' @slot constr_id (Internal) A unique integer identification number used internally.
#' @name Constraint-class
#' @aliases Constraint
#' @rdname Constraint-class
setClass("Constraint", representation(constr_id = "integer"), contains = "VIRTUAL")

setMethod("initialize", "Constraint", function(.Object, constr_id = get_id()) {
  .Object@constr_id <- constr_id
  .Object
})

setClassUnion("ListORConstr", c("list", "Constraint"))

# Helper function since syntax is different for LinOp (list) vs. Constraint object
setMethod("constr_id", "ListORConstr", function(object) {
  if(is.list(object))
    object$constr_id
  else
    object@id
})

#'
#' The BoolConstr class.
#'
#' This class represents a boolean constraint, \eqn{X_{ij} \in \{0,1\}} for all \eqn{i,j}.
#' 
#' @slot constr_id (Internal) A unique integer identification number used internally.
#' @slot lin_op A list representing the linear operator equal to the \code{.noncvx_var}.
#' @slot .noncvx_var (Internal) A list representing the variable constrained to be elementwise boolean.
#' @name BoolConstr-class
#' @aliases BoolConstr
#' @rdname BoolConstr-class
.BoolConstr <- setClass("BoolConstr", representation(lin_op = "list", .noncvx_var = "list"),
                                      prototype(.noncvx_var = NULL),
                                      validity = function(object) {
                                        if(!is.null(object@.noncvx_var))
                                          stop("[BoolConstr: .noncvx_var] .noncvx_var is an internal slot and should not be set by user")
                                        return(TRUE)
                                      }, contains = "Constraint")

#' @param lin_op A list representing the linear operator equal to the \code{.noncvx_var}.
#' @rdname BoolConstr-class
BoolConstr <- function(lin_op) { .BoolConstr(lin_op = lin_op) }

setMethod("initialize", "BoolConstr", function(.Object, ..., lin_op, .noncvx_var) {
  .Object@lin_op <- lin_op
  # Create a new non-convex variable unless the lin_op is a variable
  if(lin_op$type == VARIABLE)
    .Object@.noncvx_var <- .Object@lin_op
  else
    .Object@.noncvx_var <- create_var(.Object@lin_op$size)
  callNextMethod(.Object, ...)
})

#' @param eq_constr A list of the equality constraints in the canonical problem.
#' @param leq_constr A list of the inequality constraints in the canonical problem.
#' @param dims A list with the dimensions of the conic constraints.
#' @param solver A string representing the solver to be called.
#' @describeIn BoolConstr Format SDP constraints as inequalities for the solver.
setMethod("format_constr", "BoolConstr", function(object, eq_constr, leq_constr, dims, solver) {
  .format <- function(object) {
    eq_constr <- list()
    # If a noncvx var was created, add an equality constraint
    if(!identical(object@.noncvx_var, object@lin_op))
      eq_constr <- c(eq_constr, create_eq(object@lin_op, object@.noncvx_var))
    list(eq_constr, list())
  }

  new_eq <- .format(object)[[1]]
  # If an equality constraint was introduced, update eq_constr and dims
  if(length(new_eq) > 0) {
    eq_constr <- c(eq_constr, new_eq)
    dims[[EQ_DIM]] <- dims[[EQ_DIM]] + prod(size(object))
  }

  # Record the .noncvx_var id
  bool_id <- get_expr_vars(object@.noncvx_var)[[1]][[1]]
  constr_type <- constr_type(object)
  dims[[constr_type]] <- c(dims[[constr_type]], bool_id)
  list(eq_constr = eq_constr, leq_constr = leq_constr, dims = dims)
})

setMethod("constr_type", "BoolConstr", function(object) { BOOL_IDS })

#' @param object A \linkS4class{BoolConstr} object.
#' @describeIn BoolConstr The dimensions of the semidefinite cone.
setMethod("size", "BoolConstr", function(object) { object@lin_op$size })

#'
#' The IntConstr class.
#'
#' This class represents an integer constraint, \eqn{X_{ij} \in \mathbf{Z}} for all \eqn{i,j}.
#' 
#' @slot constr_id (Internal) A unique integer identification number used internally.
#' @slot lin_op A list representing the linear operator equal to the \code{.noncvx_var}.
#' @slot .noncvx_var (Internal) A list representing the variable constrained to be elementwise integer.
#' @name IntConstr-class
#' @aliases IntConstr
#' @rdname IntConstr-class
.IntConstr <- setClass("IntConstr", contains = "BoolConstr")

#' @param lin_op A list representing the linear operator equal to the \code{.noncvx_var}.
#' @rdname IntConstr-class
IntConstr <- function(lin_op) { .IntConstr(lin_op = lin_op) }

setMethod("constr_type", "IntConstr", function(object) { INT_IDS })

#'
#' The LeqConstraint class.
#'
#' This class represents a \eqn{\leq} inequality constraint.
#' 
#' @slot constr_id (Internal) A unique integer identification number used internally.
#' @slot lh_exp An \linkS4class{Expression}, numeric element, vector, or matrix representing the left-hand side of the inequality.
#' @slot rh_exp An \linkS4class{Expression}, numeric element, vector, or matrix representing the right-hand side of the inequality.
#' @slot args (Internal) A list that holds \code{lh_exp} and \code{rh_exp} for internal use.
#' @slot .expr (Internal) An \linkS4class{Expression} representing \code{lh_exp - rh_exp} for internal use.
#' @slot dual_variable (Internal) A \linkS4class{Variable} representing the dual variable associated with the constraint.
#' @name LeqConstraint-class
#' @aliases LeqConstraint
#' @rdname LeqConstraint-class
.LeqConstraint <- setClass("LeqConstraint", representation(lh_exp = "ConstValORExpr", rh_exp = "ConstValORExpr", args = "list", .expr = "Expression", dual_variable = "Variable"),
                           prototype(args = list(), .expr = NULL, dual_variable = Variable()),
                           validity = function(object) {
                             if(!is.null(object@.expr))
                               stop("[LeqConstraint: .expr] .expr is an internal slot and should not be set by user")
                             return(TRUE)
                           },
                            contains = c("Canonical", "Constraint"))

#' @param lh_exp An \linkS4class{Expression}, numeric element, vector, or matrix representing the left-hand side of the inequality.
#' @param rh_exp An \linkS4class{Expression}, numeric element, vector, or matrix representing the right-hand side of the inequality.
#' @rdname LeqConstraint-class
LeqConstraint <- function(lh_exp, rh_exp) { .LeqConstraint(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("initialize", "LeqConstraint", definition = function(.Object, ..., lh_exp, rh_exp, args = list(), .expr, dual_variable = Variable()) {
  .Object@lh_exp <- lh_exp
  .Object@rh_exp <- rh_exp
  .Object@args <- list(lh_exp, rh_exp)
  .Object@.expr <- lh_exp - rh_exp
  size <- size(.Object@.expr)
  .Object@dual_variable <- Variable(rows = size[1], cols = size[2])
  callNextMethod(.Object, ...)
})

setMethod("show", "LeqConstraint", function(object) {
  arg1 <- paste(as.character(object@args[[1]]), collapse = ", ")
  arg2 <- paste(as.character(object@args[[2]]), collapse = ", ")
  cat(class(object), "(", arg1, ", ", arg2, ")", sep = "")
})

#' @param x,object A \linkS4class{LeqConstraint} object.
#' @rdname LeqConstraint-class
setMethod("as.character", "LeqConstraint", function(x) {
  arg1 <- paste(as.character(x@args[[1]]), collapse = ", ")
  arg2 <- paste(as.character(x@args[[2]]), collapse = ", ")
  paste(arg1, "<=", arg2)   # TODO: Add OP_NAME parameter to LeqConstraint
})

#' @describeIn LeqConstraint The \code{constr_id} of the constraint.
setMethod("id", "LeqConstraint", function(object) { object@constr_id })

#' @describeIn LeqConstraint The size of the left-hand expression minus the right-hand expression.
setMethod("size", "LeqConstraint", function(object) { size(object@.expr) })

#' @describeIn LeqConstraint The constraint is DCP if the left-hand expression is convex and the right-hand expression is concave.
setMethod("is_dcp", "LeqConstraint", function(object) { is_convex(object@.expr) })

#' @describeIn LeqConstraint The graph implementation of the object. Marks the top level constraint as the \code{dual_holder} so the dual value will be saved to the \linkS4class{LeqConstraint}.
setMethod("canonicalize", "LeqConstraint", function(object) {
  canon <- canonical_form(object@.expr)
  dual_holder <- create_leq(canon[[1]], constr_id = object@constr_id)
  list(NA, c(canon[[2]], list(dual_holder)))
})

#' @describeIn LeqConstraint List of \linkS4class{Variable} objects in the constraint.
setMethod("variables", "LeqConstraint", function(object) { variables(object@.expr) })

#' @describeIn LeqConstraint List of \linkS4class{Parameter} objects in the constraint.
setMethod("parameters", "LeqConstraint", function(object) { parameters(object@.expr) })

#' @describeIn LeqConstraint List of \linkS4class{Constant} objects in the constraint.
setMethod("constants", "LeqConstraint", function(object) { constants(object@.expr) })

#' @describeIn LeqConstraint The elementwise maximum of the left-hand expression minus the right-hand expression, i.e. \code{max_elemwise(lh_exp - rh_exp, 0)}.
setMethod("residual", "LeqConstraint", function(object) { MaxElemwise(object@.expr, 0) })

#' @describeIn LeqConstraint A logical value indicating whether the constraint holds. Tolerance is currently set at \code{1e-4}.
setMethod("value", "LeqConstraint", function(object) {
  resid <- value(residual(object))
  if(length(resid) == 1 && is.na(resid))
    return(NA)
  else
    return(all(resid <= 1e-4))   # TODO: Add TOLERANCE parameter to LeqConstraint
})

#' @describeIn LeqConstraint A matrix representing the amount by which the constraint is off, i.e. the numeric value of the residual expression.
setMethod("violation", "LeqConstraint", function(object) { value(residual(object)) })

#' @describeIn LeqConstraint The value of the dual variable.
setMethod("dual_value", "LeqConstraint", function(object) { value(object@dual_variable) })

# Set the value of the dual variable for the constraint's parent.
setMethod("save_value", "LeqConstraint", function(object, value) { save_value(object@dual_variable, value) })

#'
#' The EqConstraint class.
#'
#' This class represents a equality constraint.
#' 
#' @slot constr_id (Internal) A unique integer identification number used internally.
#' @slot lh_exp An \linkS4class{Expression}, numeric element, vector, or matrix representing the left-hand side of the inequality.
#' @slot rh_exp An \linkS4class{Expression}, numeric element, vector, or matrix representing the right-hand side of the inequality.
#' @slot args (Internal) A list that holds \code{lh_exp} and \code{rh_exp} for internal use.
#' @slot .expr (Internal) An \linkS4class{Expression} representing \code{lh_exp - rh_exp} for internal use.
#' @slot dual_variable (Internal) A \linkS4class{Variable} representing the dual variable associated with the constraint.
#' @name EqConstraint-class
#' @aliases EqConstraint
#' @rdname EqConstraint-class
.EqConstraint <- setClass("EqConstraint", contains = "LeqConstraint")

#' @param lh_exp An \linkS4class{Expression}, numeric element, vector, or matrix representing the left-hand side of the inequality.
#' @param rh_exp An \linkS4class{Expression}, numeric element, vector, or matrix representing the right-hand side of the inequality.
#' @rdname EqConstraint-class
EqConstraint <- function(lh_exp, rh_exp) { .EqConstraint(lh_exp = lh_exp, rh_exp = rh_exp) }

#' @param object An \linkS4class{EqConstraint} object.
#' @describeIn EqConstraint The constraint is DCP if the left-hand and right-hand expressions are affine.
setMethod("is_dcp", "EqConstraint", function(object) { is_affine(object@.expr) })

#' @describeIn EqConstraint The absolute value of the left-hand minus the right-hand expression, i.e. \code{abs(lh_exp - rh_exp)}.
setMethod("residual", "EqConstraint", function(object) { abs(object@.expr) })

#' @describeIn EqConstraint The graph implementation of the object. Marks the top level constraint as the \code{dual_holder} so the dual value will be saved to the \linkS4class{EqConstraint}.
setMethod("canonicalize", "EqConstraint", function(object) {
  canon <- canonical_form(object@.expr)
  dual_holder <- create_eq(canon[[1]], constr_id = object@constr_id)
  list(NA, c(canon[[2]], list(dual_holder)))
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
.NonlinearConstraint <- setClass("NonlinearConstraint", representation(f = "function", vars_ = "list", .x_size = "numeric"),
                                 prototype(.x_size = NA_integer_), contains = "Constraint")

#' @param f A nonlinear function.
#' @param vars_ A list of variables involved in the function.
#' @rdname NonlinearConstraint-class
NonlinearConstraint <- function(f, vars_) { .NonlinearConstraint(f = f, vars_ = vars_) }

setMethod("initialize", "NonlinearConstraint", function(.Object, ..., f, vars_) {
  .Object@f <- f
  .Object@vars_ <- vars_

  # The shape of vars_ in f(vars_)
  cols <- size(.Object@vars_[[1]])[2]
  rows <- sum(sapply(.Object@vars_, function(var) { size(var)[1] }))
  .Object@.x_size <- c(rows*cols, 1)
  callNextMethod(.Object, ...)
})

#' @param object A \linkS4class{NonlinearConstraint} object.
#' @describeIn NonlinearConstraint The variables involved in the function in order, i.e. \code{f(vars_) = f(vstack(variables))}.
setMethod("variables", "NonlinearConstraint", function(object) { object@vars_ })

setMethod("block_add", "NonlinearConstraint", function(object, mat, block, vert_offset, horiz_offset, rows, cols, vert_step = 1, horiz_step = 1) {
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
  for(var in variables(object)) {
    var_size <- prod(size(var))
    var_x0 <- x0[offset:(offset + var_size)]
    big_x <- block_add(object, big_x, var_x0, var_offsets[get_data(var)], 0, var_size, 1)
    offset <- offset + var_size
  }
  big_x
})

# Place \code{Df} in the gradient of all functions.
setMethod("place_Df", "NonlinearConstraint", function(object, big_Df, Df, var_offsets, vert_offset) {
  horiz_offset <- 0
  for(var in variables(object)) {
    var_size <- prod(size(var))
    var_Df <- Df[, horiz_offset:(horiz_offset + var_size)]
    big_Df <- block_add(object, big_Df, var_Df, vert_offset, var_offsets[get_data(var)], prod(size(object)), var_size)
    horiz_offset <- horiz_offset + var_size
  }
  big_Df
})

# Place \code{H} in the Hessian of all functions.
setMethod("place_H", "NonlinearConstraint", function(object, big_H, H, var_offsets) {
  offset <- 0
  for(var in variables(object)) {
    var_size <- prod(size(var))
    var_H <- H[offset:(offset + var_size), offset:(offset + var_size)]
    big_H <- block_add(object, big_H, var_H, var_offsets[get_data(var)], var_offsets[get_data(var)], var_size, var_size)
    offset <- offset + var_size
  }
  big_H
})

# Extract the function variables from the vector \code{x} of all variables.
setMethod("extract_variables", "NonlinearConstraint", function(object, x, var_offsets) {
  local_x <- matrix(0, nrow = object@.x_size[1], ncol = object@.x_size[2])
  offset <- 0
  for(var in variables(object)) {
    var_size <- prod(size(var))
    value <- x[var_offsets[get_data(var)]:(var_offsets[get_data(var)] + var_size)]
    local_x <- block_add(object, local_x, value, offset, 0, var_size, 1)
    offset <- offset + var_size
  }
  local_x
})

# TODO: This should inherit from NonlinearConstraint once that is complete.
#'
#' The ExpCone class.
#'
#' This class represents a reformulated exponential cone constraint operating elementwise on \eqn{a, b, c}.
#' 
#' Original cone:
#' \deqn{
#' K = \{(a,b,c) | b > 0, be^{a/b} \leq c\} \cup \{(a,b,c) | a \leq 0, b = 0, c \geq 0\}
#' }
#' Reformulated cone:
#' \deqn{
#' K = \{(a,b,c) | b, c > 0, b\log(b) + a \leq b\log(c)\} \cup \{(a,b,c) | a \leq 0, b = 0, c \geq 0\}
#' }
#' 
#' @slot constr_id (Internal) A unique integer identification number used internally.
#' @slot a The variable \eqn{a} in the exponential cone.
#' @slot b The variable \eqn{b} in the exponential cone.
#' @slot c The variable \eqn{c} in the exponential cone.
#' @name ExpCone-class
#' @aliases ExpCone
#' @rdname ExpCone-class
.ExpCone <- setClass("ExpCone", representation(a = "list", b = "list", c = "list"), contains = "Constraint")

#' @param a The variable \eqn{a} in the exponential cone.
#' @param b The variable \eqn{b} in the exponential cone.
#' @param c The variable \eqn{c} in the exponential cone.
#' @rdname ExpCone-class
ExpCone <- function(a, b, c) { .ExpCone(a = a, b = b, c = c) }

setMethod("initialize", "ExpCone", function(.Object, ..., a, b, c) {
  .Object@a <- a
  .Object@b <- b
  .Object@c <- c
  callNextMethod(.Object, ...)
})

# TODO: Is this the correct size method for the exponential cone class?
#' @param x,object A \linkS4class{ExpCone} object.
#' @describeIn ExpCone The size of the \code{x} argument.
setMethod("size", "ExpCone", function(object) { size(object@a) })

#' @rdname ExpCone-class
setMethod("as.character", "ExpCone", function(x) {
  paste("ExpCone(", as.character(x@a), ", ", as.character(x@b), ", ", as.character(x@c), ")", sep = "")
})

#' @describeIn ExpCone List of \linkS4class{Variable} objects in the exponential cone.
setMethod("variables", "ExpCone", function(object) { list(object@a, object@b, object@c) })

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
    list(constraints, list(), object)
  }

  .ecos_format <- function(object) {
    list(list(), format_elemwise(list(object@a, object@c, object@b)))
  }

  .scs_format <- function(object) {
    list(list(), format_elemwise(list(object@a, object@b, object@c)))
  }

  if(is(solver, "CVXOPT"))
    stop("CVXOPT formatting has not been implemented")
    # tmp <- .cvxopt_format(object)
    # object <- tmp[[3]]   # TODO: Need to pass along updated constraint object
    # eq_constr <- c(eq_constr, tmp[[1]])
  else if(is(solver, "SCS"))
    leq_constr <- c(leq_constr, .scs_format(object)[[2]])
  else if(is(solver, "ECOS"))
    leq_constr <- c(leq_constr, .ecos_format(object)[[2]])
  else
    stop("Solver does not support exponential cone")
  # Update dims
  dims[[EXP_DIM]] <- dims[[EXP_DIM]] + prod(size(object))
  list(eq_constr = eq_constr, leq_constr = leq_constr, dims = dims)
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
.PSDConstraint <- setClass("PSDConstraint",
                           validity = function(object) {
                             lh_exp <- object@lh_exp
                             rh_exp <- object@rh_exp
                             if(size(lh_exp)[1] != size(lh_exp)[2] || size(rh_exp)[1] != size(rh_exp)[2])
                               stop("[PSDConstraint: validation] non-square matrix in positive definite constraint")
                             return(TRUE)
                           }, contains = "LeqConstraint")

#' @param lh_exp An \linkS4class{Expression}, numeric element, vector, or matrix representing the left-hand side of the inequality.
#' @param rh_exp An \linkS4class{Expression}, numeric element, vector, or matrix representing the right-hand side of the inequality.
#' @rdname PSDConstraint-class
PSDConstraint <- function(lh_exp, rh_exp) { .PSDConstraint(lh_exp = lh_exp, rh_exp = rh_exp) }

#' @param object A \linkS4class{PSDConstraint} object.
#' @describeIn PSDConstraint The constraint is DCP if the left-hand and right-hand expressions are affine.
setMethod("is_dcp", "PSDConstraint", function(object) { is_affine(object@.expr) })

#' @describeIn PSDConstraint A \linkS4class{Expression} representing the residual of the constraint.
setMethod("residual", "PSDConstraint", function(object) {
  min_eig <- LambdaMin(object@.expr + t(object@.expr))/2
  -MinElemwise(min_eig, 0)
})

#' @describeIn PSDConstraint The graph implementation of the object. Marks the top level constraint as the \code{dual_holder} so the dual value will be saved to the \linkS4class{PSDConstraint}.
setMethod("canonicalize", "PSDConstraint", function(object) {
  canon <- canonical_form(object@.expr)
  obj <- canon[[1]]
  constraints <- canon[[2]]
  half <- create_const(0.5, c(1,1))
  symm <- lo.mul_expr(half, lo.sum_expr(list(obj, lo.transpose(obj))), obj$size)
  dual_holder <- SDP(symm, enforce_sym = FALSE, constr_id = object@constr_id)
  list(NA, c(constraints, list(dual_holder)))
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
.SOC <- setClass("SOC", representation(t = "ConstValORExpr", x_elems = "list"),
                        prototype(t = NA_real_, x_elems = list()), contains = "Constraint")

#' @param t The scalar part of the second-order constraint.
#' @param x_elems A list containing the elements of the vector part of the constraint.
#' @rdname SOC-class
SOC <- function(t, x_elems) { .SOC(t = t, x_elems = x_elems) }

setMethod("initialize", "SOC", function(.Object, ..., t, x_elems) {
  .Object@t <- t
  .Object@x_elems <- x_elems
  callNextMethod(.Object, ...)
})

#' @param x,object A \linkS4class{SOC} object.
#' @rdname SOC-class
setMethod("as.character", "SOC", function(x) {
  paste("SOC(", as.character(x@t), ", <", paste(lapply(x@x_elems, function(elem) { as.character(elem) }), collapse = ", "), ">)", sep = "")
})

#' @param eq_constr A list of the equality constraints in the canonical problem.
#' @param leq_constr A list of the inequality constraints in the canonical problem.
#' @param dims A list with the dimensions of the conic constraints.
#' @param solver A string representing the solver to be called.
#' @describeIn SOC Format SOC constraints as inequalities for the solver.
setMethod("format_constr", "SOC", function(object, eq_constr, leq_constr, dims, solver) {
  .format <- function(object) {
    leq_constr <- lapply(object@x_elems, function(elem) { create_geq(elem) })
    leq_constr <- c(list(create_geq(object@t)), leq_constr)
    list(list(), leq_constr)
  }

  leq_constr <- c(leq_constr, .format(object)[[2]])
  dims[[SOC_DIM]] <- c(dims[[SOC_DIM]], size(object)[1])
  list(eq_constr = eq_constr, leq_constr = leq_constr, dims = dims)
})

#' @describeIn SOC The dimensions of the second-order cone.
setMethod("size", "SOC", function(object) {
  sizes <- sapply(object@x_elems, function(elem) { prod(size(elem)) })
  rows <- sum(sizes) + 1
  c(rows, 1)
})

#'
#' The SDP class.
#'
#' This class represents a semidefinite cone constraint, the set of all symmetric matrices such that the quadratic form \eqn{x^TAx} is non-negative for all \eqn{x}.
#' \deqn{
#' \{\mbox{symmetric } A | x^TAx \geq 0 \mbox{ for all } x\}
#' }
#' 
#' @slot constr_id (Internal) A unique integer identification number used internally.
#' @slot A The matrix variable constrained to be semidefinite.
#' @slot enforce_sym A logical value indicating whether symmetry constraints should be added.
#' @name SDP-class
#' @aliases SDP
#' @rdname SDP-class
.SDP <- setClass("SDP", representation(A = "ListORExpr", enforce_sym = "logical"),
                       prototype(enforce_sym = TRUE), contains = "Constraint")

#' @param constr_id (Internal) A unique integer identification number used internally.
#' @param A The matrix variable constrained to be semidefinite.
#' @param enforce_sym A logical value indicating whether symmetry constraints should be added.
#' @rdname SDP-class
SDP <- function(A, enforce_sym = TRUE, constr_id) {
  if(missing(constr_id))
    .SDP(A = A, enforce_sym = enforce_sym)
  else
    .SDP(A = A, enforce_sym = enforce_sym, constr_id = constr_id)
}

setMethod("initialize", "SDP", function(.Object, ..., A, enforce_sym = TRUE) {
  .Object@A <- A
  .Object@enforce_sym <- enforce_sym
  callNextMethod(.Object, ...)
})

#' @param x,object A \linkS4class{SDP} object.
#' @rdname SDP-class
setMethod("as.character", "SDP", function(x) { paste("SDP(", x@A, ")", sep = "") })

#
# Scaled Lower Triangle
# 
# Returns a linear operator representing the lower triangular entries.
# Scales the strictly lower triangular entries by \eqn{\sqrt{2}} as required by SCS.
#
# @param object A \linkS4class{SDP} object.
# @return A list representing the linear operator.
# @rdname scaled_lower_tri-int
.scaled_lower_tri <- function(object) {
  rows <- size(object)[1]
  cols <- rows
  entries <- rows*(cols + 1)/2

  val_arr <- c()
  row_arr <- c()
  col_arr <- c()
  count <- 1
  for(j in 1:cols) {
    for(i in 1:rows) {
      if(j <= i) {
        # Index in the original matrix
        col_arr <- c(col_arr, (j-1)*rows + i)

        # Index in the extracted vector
        row_arr <- c(row_arr, count)
        if(j == i)
          val_arr <- c(val_arr, 1.0)
        else
          val_arr <- c(val_arr, sqrt(2))
        count <- count + 1
      }
    }
  }

  size <- c(entries, rows*cols)
  coeff <- sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = size)
  coeff <- create_const(coeff, size, sparse = TRUE)
  vect <- lo.reshape(object@A, c(rows*cols, 1))
  lo.mul_expr(coeff, vect, c(entries, 1))
}

#
# Get Equality Constraint
# 
# Returns the equality constraints for the SDP constraint.
# 
# @param object A \linkS4class{SDP} object.
# @return A list representing the equality constraint linear operator.
# @rdname get_eq_constr-int
.get_eq_constr <- function(object) {
  upper_tri <- lo.upper_tri(object@A)
  lower_tri <- lo.upper_tri(lo.transpose(object@A))
  create_eq(upper_tri, lower_tri)
}

#' @describeIn SDP The dimensions of the semidefinite cone.
setMethod("size", "SDP", function(object) { size(object@A) })

#' @param eq_constr A list of the equality constraints in the canonical problem.
#' @param leq_constr A list of the inequality constraints in the canonical problem.
#' @param dims A list with the dimensions of the conic constraints.
#' @param solver A string representing the solver to be called.
#' @describeIn SDP Format SDP constraints as inequalities for the solver.
setMethod("format_constr", "SDP", function(object, eq_constr, leq_constr, dims, solver) {
  .scs_format <- function(object) {
    eq_constr <- .get_eq_constr(object)
    term <- .scaled_lower_tri(object)
    if(is.na(object@constr_id))
      leq_constr <- create_geq(term)
    else
      leq_constr <- create_geq(term, constr_id = object@constr_id)
    list(list(eq_constr), list(leq_constr))
  }

  if(is(solver, "CVXOPT") || is(solver, "MOSEK"))
    stop("Formatting unimplemented for CVXOPT and MOSEK")
  else if(is(solver, "SCS")) {
    scs_form <- .scs_format(object)
    new_eq_constr <- scs_form[[1]]
    new_leq_constr <- scs_form[[2]]
  } else
    stop("Solver does not support positive semidefinite cone")

  if(object@enforce_sym) {
    # upper_tri(A) == upper_tri(t(A))
    eq_constr <- c(eq_constr, new_eq_constr)

    # Update dims
    size <- size(object)
    dims[[EQ_DIM]] <- dims[[EQ_DIM]] + floor(size[1]*(size[2] - 1)/2)
  }
  # 0 <= A
  leq_constr <- c(leq_constr, new_leq_constr)

  # Update dims
  dims[[SDP_DIM]] <- c(dims[[SDP_DIM]], size(object)[1])
  list(eq_constr = eq_constr, leq_constr = leq_constr, dims = dims)
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
.SOCAxis <- setClass("SOCAxis", representation(axis = "numeric"),
                    validity = function(object) {
                      if(size(object@t)[2] != 1)
                        stop("[SOCAxis: t] t must have second dimension equal to 1")
                      if(!(length(object@axis) == 1 && object@axis %in% c(1,2)))
                        stop("[SOCAxis: axis] axis must equal 1 (row) or 2 (column)")
                      return(TRUE)
                    }, contains = "SOC")


#' @param t The scalar part of the second-order constraint.
#' @param X A matrix whose rows/columns are each a cone.
#' @param axis The dimension across which to take the slice: \code{1} indicates rows, and \code{2} indicates columns.
#' @rdname SOCAxis-class
SOCAxis <- function(t, X, axis) { .SOCAxis(t = t, x_elems = list(X), axis = axis) }

setMethod("initialize", "SOCAxis", function(.Object, ..., axis) {
  .Object@axis <- axis
  callNextMethod(.Object, ...)
})

#' @param x,object A \linkS4class{SOCAxis} object.
#' @rdname SOCAxis-class
setMethod("as.character", "SOCAxis", function(x) {
  paste("SOCAxis(", as.character(x@t), ", ", as.character(x@x_elems[[1]]), ", <", paste(x@axis, collapse = ", "), ">)", sep = "")
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
setMethod("num_cones", "SOCAxis", function(object) { size(object@t)[1] })

#' @describeIn SOCAxis The dimensions of a single cone.
setMethod("cone_size", "SOCAxis", function(object) {
  if(object@axis == 1)   # Return ncols if applying along each row
    c(1 + size(object@x_elems[[1]])[2], 1)
  else if(object@axis == 2)   # Return nrows if applying along each column
    c(1 + size(object@x_elems[[1]])[1], 1)
  else
    stop("Invalid axis ", object@axis)
})

#' @describeIn SOCAxis The dimensions of the (elementwise) second-order cones.
setMethod("size", "SOCAxis", function(object) {
  cone_size <- cone_size(object)
  lapply(1:num_cones(object), function(i) { cone_size })
})
