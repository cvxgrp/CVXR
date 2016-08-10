##setClass("Constraint", representation(id = "character"), prototype(id = uuid::UUIDgenerate()), contains = "VIRTUAL")
setClass("Constraint", representation(constr_id = "integer"), prototype(constr_id = get_id()), contains = "VIRTUAL")
setMethod("id", "Constraint", function(object) { object@constr_id })

.BoolConstr <- setClass("BoolConstr", representation(lin_op = "list", .noncvx_var = "list"),
                                      prototype(.noncvx_var = NULL),
                                      validity = function(object) {
                                        if(!is.null(object@.noncvx_var))
                                          stop("[BoolConstr: .noncvx_var] .noncvx_var is an internal slot and should not be set by user")
                                        return(TRUE)
                                      }, contains = "Constraint")
BoolConstr <- function(lin_op) { .BoolConstr(lin_op = lin_op) }

setMethod("initialize", "BoolConstr", function(.Object, ..., lin_op, .noncvx_var) {
  .Object@lin_op <- lin_op
  if(lin_op$type == VARIABLE)
    .Object@.noncvx_var <- .Object@lin_op
  else
    .Object@.noncvx_var <- create_var(.Object@lin_op$size)
  callNextMethod(.Object, ...)
})

setMethod("format_constr", "BoolConstr", function(object, eq_constr, leq_constr, dims, solver) {
  .format <- function(object) {
    eq_constr <- list()
    if(object@noncvx_var != object@lin_op)
      eq_constr <- c(eq_constr, create_eq(object@lin_op, object@noncvx_var))
    list(eq_constr, list())
  }

  new_eq <- .format(object)[[1]]
  if(length(new_eq) > 0) {
    eq_constr <- c(eq_constr, new_eq)
    dims[EQ_DIM] <- c(dims[EQ_DIM], size(object)[1] * size(object)[2])
  }

  bool_id <- get_expr_vars(object@noncvx_var)[[1]][[1]]
  CONSTR_TYPE <- constr_type(object)
  dims[CONSTR_TYPE] <- c(dims[CONSTR_TYPE], bool_id)
})

setMethod("constr_type", "BoolConstr", function(object) { BOOL_IDS })
setMethod("size", "BoolConstr", function(object) { .Object@lin_op$size })

IntConstr <- setClass("IntConstr", contains = "BoolConstr")

setMethod("constr_type", "IntConstr", function(object) { INT_IDS })

.LeqConstraint <- setClass("LeqConstraint", representation(lh_exp = "Expression", rh_exp = "Expression", .expr = "Expression", dual_variable = "Variable"),
                           prototype(.expr = NULL, dual_variable = new("Variable")),
                           validity = function(object) {
                             if(!is.null(object@.expr))
                               stop("[LeqConstraint: .expr] .expr is an internal slot and should not be set by user")
                             return(TRUE)
                           },
                            contains = "Constraint")
LeqConstraint <- function(lh_exp, rh_exp) { .LeqConstraint(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("initialize", "LeqConstraint", definition = function(.Object, ..., lh_exp, rh_exp, .expr) {
  .Object@lh_exp <- lh_exp
  .Object@rh_exp <- rh_exp
  .Object@.expr <- lh_exp - rh_exp
  .Object@dual_variable <- Variable(size(.Object@.expr))
  callNextMethod(.Object, ...)
})

setMethod("size", "LeqConstraint", function(object) { size(object@.expr) })
setMethod("is_dcp", "LeqConstraint", function(object) { is_convex(object@.expr) })
setMethod("canonical_form", "LeqConstraint", function(object) { canonicalize(object) })

setMethod("canonicalize", "LeqConstraint", function(object) {
  canon <- canonical_form(object@.expr)
  dual <- create_leq(canon[[1]], constr_id = object@constr_id)
  list(NA, c(canon[[2]], list(dual)))
})

setMethod("variables", "LeqConstraint", function(object) { variables(object@.expr) })
setMethod("dual_variable", "LeqConstraint", function(object) { value(object@dual_variable) })
setMethod("parameters", "LeqConstraint", function(object) { parameters(object@.expr) })
setMethod("value", "LeqConstraint", function(object) {
  if(is.na(value(object@.expr)))
    NA
  else
    all(value(object@.expr) <= 1e-6)   # TODO: Add tolerance parameter to LeqConstraint
})

setMethod("violation", "LeqConstraint", function(object) {
  if(is.na(value(object@.expr)))
    NA
  else
    max(value(object@.expr), 0)
})

NonlinearConstraint <- setClass("NonlinearConstraint", representation(f = "Expression", vars = "list", .x_size = "numeric"),
                                prototype(vars = list(), .x_size = NA_real_),
                                validity = function(object) {
                                  if(!is.na(.x_size))
                                    stop("[NonlinearConstraint: .x_size] .x_size is an internal slot and should not be set by user")
                                  is_var <- sapply(vars, function(var) { is(var, "Variable") })
                                  if(!all(is_var))
                                    stop("[NonlinearConstraint: vars] vars must be a list of Variable objects")
                                  return(TRUE)
                                  }, contains = "Constraint")

setMethod("initialize", "NonlinearConstraint", function(.Object, ..., f, vars) {
  object@f <- f
  object@vars <- vars
  cols <- size(object@vars[[1]])[2]
  rows <- sum(sapply(object@vars, function(var) { size(var)[1] }))
  object@.x_size <- c(rows*cols, 1)
  callNextMethod(.Object, ...)
})

setMethod("variables", "NonlinearConstraint", function(object) { object@vars })

EqConstraint <- setClass("EqConstraint", contains = "LeqConstraint")

setMethod("is_dcp", "EqConstraint", function(object) { is_affine(object@.expr) })
setMethod("value", "EqConstraint", function(object) {
  if(is.na(value(object@.expr)))
    NA
  else
    all(abs(value(object@.expr)) <= 1e-6)   # TODO: Set tolerance as parameter of class
})

setMethod("violation", "EqConstraint", function(object) {
  if(is.na(value(object@.expr)))
    NA
  else
    abs(value(object@.expr))
})

setMethod("canonicalize", "EqConstraint", function(object) {
  canon <- canonical_form(object@.expr)
  dual <- create_eq(canon[[1]], constr_id = object@constr_id)
  c(NA, c(canon[[2]], dual))
})

.ExpCone <- setClass("ExpCone", representation(x = "ConstValORExpr", y = "ConstValORExpr", z = "ConstValORExpr"), contains = "NonlinearConstraint")
ExpCone <- function(x, y, z) { .ExpCone(x = x, y = y, z = z) }

setMethod("initialize", "ExpCone", function(.Object, ..., x, y, z) {
  .Object@x <- x
  .Object@y <- y
  .Object@z <- z
  callNextMethod("ExpCone", .Object, ..., vars = list(.Object@x, .Object@y, .Object@z))
})

setMethod("size", "ExpCone", function(object) { size(object@x) })

.PSDConstraint <- setClass("PSDConstraint", contains = "LeqConstraint",
                           validity = function(object) {
                             lh_exp <- object@lh_exp
                             rh_exp <- object@rh_exp
                             if(size(lh_exp)[1] != size(lh_exp)[2] || size(rh_exp)[1] != size(rh_exp)[2])
                               stop("[PSDConstraint: validation] non-square matrix in positive definite constraint")
                             return(TRUE)
                           })
PSDConstraint <- function(lh_exp, rh_exp) { .PSDConstraint(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("is_dcp", "PSDConstraint", function(object) { is_affine(object@.expr) })
setMethod("value", "PSDConstraint", function(object) {
  if(is.na(value(object@.expr)))
    NA
  else {
    mat <- value(object@.expr)
    w <- eigen(mat + t(mat))
    min(w$values)/2 >= -1e-6
  }
})

setMethod("violation", "PSDConstraint", function(object) {
  if(is.na(value(object@.expr)))
    NA
  else {
    mat <- value(object@.expr)
    w <- eigen(mat + t(mat))
    -min(min(w$values)/2, 0)
  }
})

setMethod("canonicalize", "PSDConstraint", function(object) {
  canon <- canonical_form(object@.expr)
  obj <- canon[[1]]
  constraints <- canon[[2]]
  half <- create_const(0.5, c(1,1))
  symm <- mul_expr(half, sum_expr(list(obj, transpose(obj))), size(obj))
  dual_holder <- SDP(symm, enforce_sym = FALSE, constr_id = object@constr_id)
  list(NA, c(constraints, dual_holder))
})

.SOC <- setClass("SOC", representation(t = "ConstValORExpr", x_elems = "list"),
                       prototype(t = NA_real_, x_elems = list()), contains = "Constraint")
SOC <- function(t, x_elems) { .SOC(t = t, x_elems = x_elems) }

setMethod("show", "SOC", function(object) { cat("SOC(", as.character(object@t), ", <", paste(lapply(object@x_elems, function(x) { as.character(x) }), collapse = ", "), ">)", sep = "") })

setMethod("format_constr", "SOC", function(object) {
  .format <- function(object) {
    leq_constr <- lapply(object@x_elems, function(elem) { create_geq(elem) })
    leq_constr <- c(create_geq(object@t), leq_constr)
    list(list(), leq_constr)
  }

  leq_constr <- c(leq_constr, .format(object)[[2]])
  dims[SOC_DIM] <- c(dims[SOC_DIM], size(object)[1])
})

setMethod("size", "SOC", function(object) {
  sizes <- sapply(object@x_elems, function(elem) {
    if(is.list(elem))
      prod(elem$size)
    else
      prod(size(elem))
  })
  rows <- sum(sizes) + 1
  c(rows, 1)
})

.SDP <- setClass("SDP", representation(A = "numeric", enforce_sym = "logical"),
                       prototype(A = NA_real_, enforce_sym = TRUE), contains = "Constraint")
SDP <- function(A, enforce_sym = TRUE, constr_id) {
  if(missing(constr_id))
    .SDP(A = A, enforce_sym = enforce_sym)
  else
    .SDP(A = A, enforce_sym = enforce_sym, constr_id = constr_id)
}

setMethod("show", "SDP", function(object) { cat("SDP(", object@A, ")", sep = "") })
setMethod("size", "SDP", function(object) { size(object@A) })

SOCElemwise <- setClass("SOCElemwise", contains = "SOC")

setMethod("format_constr", "SOCElemwise", function(object, eq_constr, leq_constr, dims, solver) {
 .format <- function(object) {
   list(list(), format_elemwise(list(object@t, object@x_elems)))
 }

 leq_constr <- c(leq_constr, .format(object)[[2]])
 for(cone_size in size(object))
   dims[SOC_DIM] <- c(dims[SOC_DIM], cone_size[1])
})

setMethod("num_cones", "SOCElemwise", function(object) {
  size(object@t)[1] * size(object@t)[2]
})

setMethod("cone_size", "SOCElemwise", function(object) {
  c(1 + length(object@x_elems), 1)
})

setMethod("size", "SOCElemwise", function(object) {
  cone_size <- cone_size(object)
  lapply(1:num_cones(object), function(i) { cone_size })
})
