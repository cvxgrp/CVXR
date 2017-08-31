setClass("Constraint", representation(constr_id = "integer"), contains = "VIRTUAL")

setMethod("initialize", "Constraint", function(.Object, constr_id = get_id()) {
  .Object@constr_id <- constr_id
  .Object
})

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
  # Create a new non-convex variable unless the lin_op is a variable
  if(lin_op$type == VARIABLE)
    .Object@.noncvx_var <- .Object@lin_op
  else
    .Object@.noncvx_var <- create_var(.Object@lin_op$size)
  callNextMethod(.Object, ...)
})

setMethod("format_constr", "BoolConstr", function(object, eq_constr, leq_constr, dims, solver) {
  .format <- function(object) {
    eq_constr <- list()
    if(object@.noncvx_var != object@lin_op)
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
  dims[[BOOL_IDS]] <- c(dims[[BOOL_IDS]], bool_id)
  list(eq_constr = eq_constr, leq_constr = leq_constr, dims = dims)
})

setMethod("constr_type", "BoolConstr", function(object) { BOOL_IDS })
setMethod("size", "BoolConstr", function(object) { object@lin_op$size })

.IntConstr <- setClass("IntConstr", contains = "BoolConstr")
IntConstr <- function(lin_op) { .IntConstr(lin_op = lin_op) }
setMethod("constr_type", "IntConstr", function(object) { INT_IDS })

.LeqConstraint <- setClass("LeqConstraint", representation(lh_exp = "ConstValORExpr", rh_exp = "ConstValORExpr", args = "list", .expr = "Expression", dual_variable = "Variable"),
                           prototype(args = list(), .expr = NULL, dual_variable = Variable()),
                           validity = function(object) {
                             if(!is.null(object@.expr))
                               stop("[LeqConstraint: .expr] .expr is an internal slot and should not be set by user")
                             return(TRUE)
                           },
                            contains = c("Canonical", "Constraint"))
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
  cat(class(object), "(", as.character(object@args[[1]]), ", ", as.character(object@args[[2]]), ")", sep = "")
})

setMethod("as.character", "LeqConstraint", function(x) {
  paste(as.character(x@args[[1]]), "<=", as.character(x@args[[2]]))   # TODO: Add OP_NAME parameter to LeqConstraint
})

setMethod("id", "LeqConstraint", function(object) { object@constr_id })
setMethod("size", "LeqConstraint", function(object) { size(object@.expr) })
setMethod("is_dcp", "LeqConstraint", function(object) { is_convex(object@.expr) })
setMethod("canonicalize", "LeqConstraint", function(object) {
  canon <- canonical_form(object@.expr)
  dual_holder <- create_leq(canon[[1]], constr_id = object@constr_id)
  list(NA, c(canon[[2]], list(dual_holder)))
})

setMethod("variables", "LeqConstraint", function(object) { variables(object@.expr) })
setMethod("parameters", "LeqConstraint", function(object) { parameters(object@.expr) })
setMethod("constants", "LeqConstraint", function(object) { constants(object@.expr) })
setMethod("residual", "LeqConstraint", function(object) { MaxElemwise(object@.expr, 0) })
setMethod("value", "LeqConstraint", function(object) {
  resid <- value(residual(object))
  if(length(resid) == 1 && is.na(resid))
    return(NA)
  else
    return(all(resid <= 1e-4))   # TODO: Add TOLERANCE parameter to LeqConstraint
})
setMethod("violation", "LeqConstraint", function(object) { value(residual(object)) })
setMethod("dual_value", "LeqConstraint", function(object) { value(object@dual_variable) })
setMethod("save_value", "LeqConstraint", function(object, value) { save_value(object@dual_variable, value) })

.EqConstraint <- setClass("EqConstraint", contains = "LeqConstraint")
EqConstraint <- function(lh_exp, rh_exp) { .EqConstraint(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("is_dcp", "EqConstraint", function(object) { is_affine(object@.expr) })
setMethod("residual", "EqConstraint", function(object) { abs(object@.expr) })
setMethod("canonicalize", "EqConstraint", function(object) {
  canon <- canonical_form(object@.expr)
  dual_holder <- create_eq(canon[[1]], constr_id = object@constr_id)
  list(NA, c(canon[[2]], list(dual_holder)))
})

# TODO: Do I need the NonlinearConstraint class?
.NonlinearConstraint <- setClass("NonlinearConstraint", representation(f = "function", vars_ = "list", .x_size = "numeric"), 
                                 prototype(.x_size = NA_integer_), contains = "Constraint")
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

setMethod("variables", "NonlinearConstraint", function(object) { object@vars_ })
setMethod("block_add", "NonlinearConstraint", function(object, mat, block, vert_offset, horiz_offset, rows, cols, vert_step = 1, horiz_step = 1) {
  row_seq <- seq(vert_offset, rows + vert_offset, vert_step)
  col_seq <- seq(horiz_offset, cols + horiz_offset, horiz_step)
  mat[row_seq, col_seq] <- mat[row_seq, col_seq] + block
  mat
})

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

# TODO: This should inherit from NonlinearConstraint once that is complete
.ExpCone <- setClass("ExpCone", representation(x = "list", y = "list", z = "list"), contains = "Constraint")
ExpCone <- function(x, y, z) { .ExpCone(x = x, y = y, z = z) }

setMethod("initialize", "ExpCone", function(.Object, ..., x, y, z) {
  .Object@x <- x
  .Object@y <- y
  .Object@z <- z
  callNextMethod(.Object, ...)
})

# TODO: Is this the correct size method for the exponential cone class?
setMethod("size", "ExpCone", function(object) { size(object@x) })

setMethod("as.character", "ExpCone", function(x) {
  paste("ExpCone(", as.character(x@x), ", ", as.character(x@y), ", ", as.character(x@z), ")", sep = "")
})
setMethod("variables", "ExpCone", function(object) { list(object@x, object@y, object@z) })
setMethod("format_constr", "ExpCone", function(object, eq_constr, leq_constr, dims, solver) {
  .cvxopt_format <- function(object) {
    constraints <- list()
    for(i in 1:length(object@vars_)) {
      var <- vars_[[i]]
      if(var$type != VARIABLE) {
        lone_var <- create_var(var$size)
        constraints <- c(constraints, list(create_eq(lone_var, var)))
        object@vars_[[i]] <- lone_var
      }
    }
    list(constraints, list(), object)
  }
  
  .ecos_format <- function(object) {
    list(list(), format_elemwise(list(object@x, object@z, object@y)))
  }
  
  .scs_format <- function(object) {
    list(list(), format_elemwise(list(object@x, object@y, object@z)))
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

.PSDConstraint <- setClass("PSDConstraint",
                           validity = function(object) {
                             lh_exp <- object@lh_exp
                             rh_exp <- object@rh_exp
                             if(size(lh_exp)[1] != size(lh_exp)[2] || size(rh_exp)[1] != size(rh_exp)[2])
                               stop("[PSDConstraint: validation] non-square matrix in positive definite constraint")
                             return(TRUE)
                           }, contains = "LeqConstraint")
PSDConstraint <- function(lh_exp, rh_exp) { .PSDConstraint(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("is_dcp", "PSDConstraint", function(object) { is_affine(object@.expr) })
setMethod("residual", "PSDConstraint", function(object) {
  min_eig <- LambdaMin(object@.expr + t(object@.expr))/2
  -MinElemwise(min_eig, 0)
})

setMethod("canonicalize", "PSDConstraint", function(object) {
  canon <- canonical_form(object@.expr)
  obj <- canon[[1]]
  constraints <- canon[[2]]
  half <- create_const(0.5, c(1,1))
  symm <- mul_expr(half, sum_expr(list(obj, transpose(obj))), obj$size)
  dual_holder <- SDP(symm, enforce_sym = FALSE, constr_id = object@constr_id)
  list(NA, c(constraints, list(dual_holder)))
})

.SOC <- setClass("SOC", representation(t = "ConstValORExpr", x_elems = "list"),
                        prototype(t = NA_real_, x_elems = list()), contains = "Constraint")
SOC <- function(t, x_elems) { .SOC(t = t, x_elems = x_elems) }

setMethod("initialize", "SOC", function(.Object, ..., t, x_elems) {
  .Object@t <- t
  .Object@x_elems <- x_elems
  callNextMethod(.Object, ...)
})

setMethod("as.character", "SOC", function(x) { 
  paste("SOC(", as.character(object@t), ", <", paste(lapply(object@x_elems, function(x) { as.character(x) }), collapse = ", "), ">)", sep = "") 
})

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

setMethod("size", "SOC", function(object) {
  sizes <- sapply(object@x_elems, function(elem) { prod(size(elem)) })
  rows <- sum(sizes) + 1
  c(rows, 1)
})

.SDP <- setClass("SDP", representation(A = "ListORExpr", enforce_sym = "logical"),
                       prototype(enforce_sym = TRUE), contains = "Constraint")
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

setMethod("as.character", "SDP", function(x) { paste("SDP(", x@A, ")", sep = "") })

.scaled_lower_tri <- function(object) {
  require(Matrix)
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
  vect <- reshape(object@A, c(rows*cols, 1))
  mul_expr(coeff, vect, c(entries, 1))
}

.get_eq_constr <- function(object) {
  upper_tri <- upper_tri(object@A)
  lower_tri <- upper_tri(transpose(object@A))
  create_eq(upper_tri, lower_tri)
}

setMethod("size", "SDP", function(object) { size(object@A) })
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

# Assumes t is a vector the same length as X's rows (columns) for axis == 1 (2)
.SOCAxis <- setClass("SOCAxis", representation(axis = "numeric"), 
                    validity = function(object) {
                      if(size(object@t)[2] != 1)
                        stop("[SOCAxis: t] t must have second dimension equal to 1")
                      if(!(length(object@axis) == 1 && object@axis %in% c(1,2)))
                        stop("[SOCAxis: axis] axis must equal 1 (row) or 2 (column)")
                      return(TRUE)
                    }, contains = "SOC")
SOCAxis <- function(t, X, axis) { .SOCAxis(t = t, x_elems = list(X), axis = axis) }

setMethod("initialize", "SOCAxis", function(.Object, ..., axis) {
  .Object@axis <- axis
  callNextMethod(.Object, ...)
})

setMethod("as.character", "SOCAxis", function(x) { 
  paste("SOCAxis(", as.character(x@t), ", ", as.character(x@x_elems[[1]]), ", <", paste(x@axis, collapse = ", "), ">)", sep = "")
})

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

setMethod("num_cones", "SOCAxis", function(object) { size(object@t)[1] })
setMethod("cone_size", "SOCAxis", function(object) {
  if(object@axis == 1)   # Return ncols if applying along each row
    c(1 + size(object@x_elems[[1]])[2], 1)
  else if(object@axis == 2)   # Return nrows if applying along each column
    c(1 + size(object@x_elems[[1]])[1], 1)
  else
    stop("Invalid axis ", object@axis)
})
setMethod("size", "SOCAxis", function(object) {
  cone_size <- cone_size(object)
  lapply(1:num_cones(object), function(i) { cone_size })
})
