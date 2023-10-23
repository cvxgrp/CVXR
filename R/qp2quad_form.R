#'
#' The Qp2SymbolicQp class.
#'
#' This class reduces a quadratic problem to a problem that consists of affine
#' expressions and symbolic quadratic forms.
#'
#' @rdname Qp2SymbolicQp-class
.Qp2SymbolicQp <- setClass("Qp2SymbolicQp", contains = "Canonicalization")
Qp2SymbolicQp <- function(problem = NULL) { .Qp2SymbolicQp(problem = problem) }

setMethod("initialize", "Qp2SymbolicQp", function(.Object, ...) {
  callNextMethod(.Object, ..., canon_methods = Qp2QuadForm.CANON_METHODS)
})

Qp2SymbolicQp.accepts <- function(problem) {
  is_qpwa(expr(problem@objective)) &&
  length(intersect(c("PSD", "NSD"), convex_attributes(variables(problem)))) == 0 &&
  all(sapply(problem@constraints, function(c) {
        (inherits(c, c("NonPosConstraint", "NonNegConstraint", "IneqConstraint")) && is_pwl(expr(c))) ||
        (inherits(c, c("ZeroConstraint", "EqConstraint")) && are_args_affine(list(c)))
  }))
}

# Problems with quadratic, piecewise affine objectives, piecewise-linear constraints, inequality constraints,
# and affine equality constraints are accepted.
setMethod("accepts", signature(object = "Qp2SymbolicQp", problem = "Problem"), function(object, problem) {
    Qp2SymbolicQp.accepts(problem)
})

# Converts a QP to an even more symbolic form.
setMethod("perform", signature(object = "Qp2SymbolicQp", problem = "Problem"), function(object, problem) {
  if(!accepts(object, problem))
    stop("Cannot reduce problem to symbolic QP")
  callNextMethod(object, problem)
})

#'
#'
#' The QpMatrixStuffingConeDims class.
#'
#' This class contains a summary of cone dimensions present in constraints.
#'
#' @slot zero The dimension of the zero cone.
#' @slot nonpos The dimension of the nonpositive cone.
#' @slot exp The number of 3-dimensional exponential cones.
#' @slot soc A vector of the second-order cone dimensions.
#' @slot psd A vector the positive semidefinite cone dimensions, where the dimension of the PSD cone of k by k matrices is k.
#' @rdname QpMatrixStuffingConeDims-class
.QpMatrixStuffingConeDims <- setClass("QpMatrixStuffingConeDims", representation(constr_map = "list", zero = "numeric", nonneg = "numeric", exp = "numeric", soc = "numeric", psd = "numeric"),
                                  prototype(zero = 0, nonpos = 0, exp = 0, soc = 0, psd = 0))
QpMatrixStuffingConeDims <- function(constr_map) { .QpMatrixStuffingConeDims(constr_map = constr_map) }

setMethod("initialize", "QpMatrixStuffingConeDims", function(.Object, constr_map, zero = 0, nonpos = 0, exp = 0, soc = 0, psd = 0) {
  .Object@constr_map <- constr_map
  .Object@zero <- as.integer(sum(sapply(constr_map$Zero, size)))
  .Object@nonneg <- as.integer(sum(sapply(constr_map$NonNeg, size)))
  .Object@exp <- as.integer(sum(sapply(constr_map$ExpCone, num_cones)))
  # .Object@soc <- c()
  # for(c in constr_map$SOC)
  #   .Object@soc <- c(.Object@soc, cone_sizes(c))
  # .Object@soc <- as.integer(.Object@soc)
  .Object@soc <- as.integer(Reduce(base::c, lapply(constr_map$SOC, cone_sizes)))
  .Object@psd <- as.integer(sapply(constr_map$PSD, nrow))
  .Object
})

setMethod("show", "QpMatrixStuffingConeDims", function(object) {
  print(paste("(zero: ", x@zero, ", nonpos: ", x@nonpos, ", exp: ", x@exp, ", soc: ", x@soc, ", psd: ", x@psd, ")", sep = ""))
})

setMethod("as.character", "QpMatrixStuffingConeDims", function(x) {
  paste(x@zero, " equalities, ", x@nonpos, " inequalities, ", x@exp, " exponential cones, \nSOC constraints: ", x@soc, ", PSD constraints: ", x@psd, sep = "")
})

setMethod("$", "QpMatrixStuffingConeDims", function(x, name) {
  if(name == EQ_DIM)
    return(x@zero)
  else if(name == LEQ_DIM)
    return(x@nonpos)
  else if(name == EXP_DIM)
    return(x@exp)
  else if(name == SOC_DIM)
    return(x@soc)
  else if(name == PSD_DIM)
    return(x@psd)
  else
    stop("Unknown key ", name)
})

#'
#' The ParamQuadProg class.
#'
#' This class represents a parametrized quadratic program.
#'
#' minimize   x^T*P*x + q^T*x + d
#' subject to (in)equality_constr(A_1*x + b_1, ...)
#'            ...
#'            (in)equality_constr(A_i*x + b_i, ...)
#'
#' The constant offsets d and b are the last column of c and A, respectively.
#'
#' @rdname ParamQuadProg-class
.ParamQuadProg <- setClass("ParamQuadProg", representation(P = "numeric", q = "numeric", x = "Expression", A = "numeric", variables = "list", var_id_to_col = "list", constraints = "list", parameters = "list", param_id_to_col = "list", formatted = "logical", reduced_P = "ReducedMatORNULL", reduced_A = "ReducedMatORNULL", constr_size = "numeric", id_to_param = "list", param_id_to_size = "list", total_param_size = "numeric", id_to_var = "list"),
                           prototype(formatted = FALSE, reduced_P = NULL, reduced_A = NULL, id_to_param = list(), param_id_to_size = list(), total_param_size = NA_integer_, id_to_var = list()),
                           contains = "ParamProb")
ParamQuadProg <- function(P, q, x, A, variables, var_id_to_col, constraints, parameters, param_id_to_col, formatted = FALSE) {
  .ParamQuadProg(P = P, q = q, x = x, A = A, variables = variables, var_id_to_col = var_id_to_col, constraints = constraints, parameters = parameters, param_id_to_col = param_id_to_col, formatted = formatted)
}

setMethod("initialize", "ParamQuadProg", function(.Object, ..., P, q, x, A, variables, var_id_to_col, constraints, parameters, param_id_to_col, formatted = FALSE, reduced_P = NULL, reduced_A = NULL, id_to_param = list(), param_id_to_size = list(), total_param_size = NA_integer_, id_to_var = list()) {
  .Object@P <- P
  .Object@q <- q
  .Object@x <- x
  .Object@A <- A

  # Form a reduced representation of A and P for faster application of parameters.
  .Object@reduced_A <- ReducedMat(.Object@A, size(.Object@x))
  .Object@reduced_P <- ReducedMat(.Object@P, size(.Object@x), quad_form = TRUE)

  .Object@constraints <- constraints
  .Object@constr_size <- sum(sapply(constraints, size))
  .Object@parameters <- parameters
  .Object@param_id_to_col <- param_id_to_col
  .Object@id_to_param <- list()
  for(p in .Object@parameters)
    .Object@id_to_param[[as.character(p@id)]] <- size(p)
  .Object@total_param_size <- sum(sapply(.Object@parameters, size))

  # TODO: Technically part of inverse data.
  .Object@variables <- variables
  .Object@var_id_to_col <- var_id_to_col
  .Object@id_to_var <- list()
  for(v in .Object@variables)
    .Object@id_to_var[[as.character(v@id)]] <- v
  .Object@formatted <- formatted   # Has this param cone prog been formatted for a solver?

  callNextMethod(.Object, ...)
})

#' @param object A \linkS4class{ParamQuadProg} object.
#' @describeIn ParamQuadProg Is the problem mixed-integer?
setMethod("is_mixed_integer", "ParamQuadProg", function(object) {
  object@x@attributes$boolean || object@x@attributes$integer
})

#' @param id_to_param_value (Optional) List mapping paramter ids to values.
#' @param zero_offset (Optional) If TRUE, zero out the constant offset in the parameter vector.
#' @param keep_zeros (Optional) If TRUE, store explicit zeros in A where parameters are affected.
#' @describeIn ParamQuadProg Returns A, b after applying parameters (and reshaping).
setMethod("apply_parameters", "ParamQuadProg", function(object, id_to_param_value = NULL, zero_offset = FALSE, keep_zeros = FALSE) {
  param_value <- function(idx) {
    if(is.null(id_to_param_value) || length(id_to_param_value) == 0)
      return(value(id_to_param[[idx]]))
    else
      return(id_to_param_value[[idx]])
  }
  param_vec <- canonInterface.get_parameter_vector(object@total_param_size, object@param_id_to_col, object@param_id_to_size, param_value, zero_offset = zero_offset)

  object@reduced_P <- cache(object@reduced_P, keep_zeros)
  P <- get_matrix_from_tensor(object@reduced_P, param_vec, with_offset = FALSE)[[1]]

  qd <- canonInterface.get_matrix_from_tensor(object@q, param_vec, size(object@x), with_offset = TRUE)
  q <- qd[[1]]
  d <- qd[[2]]
  q <- as.matrix(as.vector(q))

  object@reduced_A <- cache(object@reduced_A, keep_zeros)
  Ab <- get_matrix_from_tensor(object@reduced_A, param_vec, with_offset = TRUE)
  A <- Ab[[1]]
  b <- Ab[[2]]

  return(list(P = P, q = q, d = d, AF = A, b = as.matrix(b)))
})

#' @describeIn ParamQuadProg Multiplies by Jacobian of parameter mapping.
setMethod("apply_param_jac", "ParamQuadProg", function(object, delP, delq, delA, delb, active_params = NULL) {
  stop("Unimplemented")
})

#' @describeIn ParamQuadProg Splits the solution into individual variables.
setMethod("split_solution", "ParamQuadProg", function(object, sltn, active_vars = NULL) {
  stop("Unimplemented")
})

#' @describeIn ParamQuadProg Adjoint of split_solution.
setMethod("split_solution", "ParamQuadProg", function(object, del_vars = NULL) {
  stop("Unimplemented")
})

#'
#' The QpMatrixStuffing class.
#'
#' This class fills in numeric values for the problem instance and
#' outputs a DCP-compliant minimization problem with an objective
#' of the form
#'
#' QuadForm(x, p) + t(q) %*% x
#'
#' and Zero/NonPos constraints, both of which exclusively carry
#' affine arguments
#'
#' @rdname QpMatrixStuffing-class
.QpMatrixStuffing <- setClass("QpMatrixStuffing", representation(canon_backend = "character"), prototype(canon_backend = NA_character_), contains = "MatrixStuffing")
QpMatrixStuffing <- function(canon_backend = NA_character_) { .QpMatrixStuffing(canon_backend = canon_backend) }

setMethod("initialize", "QpMatrixStuffing", function(.Object, ..., canon_backend = NA_character_) {
  .Object@canon_backend <- canon_backend
  callNextMethod(.Object, ...)
})

setMethod("accepts", signature(object = "QpMatrixStuffing", problem = "Problem"), function(object, problem) {
    inherits(problem@objective, "Minimize") &&
    is_quadratic(problem@objective) &&
    is_dcp(problem) &&
    length(convex_attributes(variables(problem))) == 0 &&
    all(sapply(problem@constraints, inherits, what = c("ZeroConstraint", "NonPosConstraint", "EqConstraint", "IneqConstraint") )) &&
    are_args_affine(problem@constraints) &&
    is_dpp(problem)
})

setMethod("stuffed_objective", signature(object = "QpMatrixStuffing", problem = "Problem", extractor = "CoeffExtractor"), function(object, problem, extractor) {
  # Extract to 0.5 * t(x) %*% P %*% x + t(q) %*% x + r
  expr <- copy(expr(problem@objective))   # TODO: Do we need to copy objective?
  Pq_params <- quad_form(extractor, expr)
  params_to_P <- Pq_params[[1]]
  params_to_q <- Pq_params[[2]]
  # Handle 0.5 factor.
  params_to_P <- 2*params_to_P

  # Concatenate all variables in one vector.
  boolint <- extract_mip_idx(variables(problem))
  boolean <- boolint[[1]]
  integer <- boolint[[2]]
  x <- Variable(extractor@x_length, boolean = boolean, integer = integer)
  return(list(params_to_P, params_to_q, x))
})

setMethod("perform", signature(object = "QpMatrixStuffing", problem = "Problem"), function(object, problem) {
  # See docs for MatrixStuffing perform function.
  inverse_data <- InverseData(problem)
  # Form the constraints
  extractor <- CoeffExtractor(inverse_data, object@canon_backend)
  tmp <- stuffed_objective(object, problem, extractor)
  params_to_P <- tmp[[1]]
  params_to_q <- tmp[[2]]
  flattened_variable <- tmp[[3]]

  # Lower equality and inequality to ZeroConstraint and NonPosConstraint.
  cons <- list()
  for(con in problem@constraints) {
    if(is(con, "EqConstraint"))
      con <- lower_equality(con)
    else if(is(con, "IneqConstraint"))
      con <- lower_ineq_to_nonneg(con)
    else if(is(con, "NonPosConstraint"))
      con <- nonpos2nonneg(con)
    cons <- c(cons, con)
  }

  # Reorder constraints to ZeroConstraint, NonNegConstraint.
  constr_map <- group_constraints(cons)
  ordered_cons <- c(constr_map$Zero, constr_map$NonNeg)
  inverse_data@cons_id_map <- list()
  for(con in ordered_cons) {
    con_id <- id(con)
    inverse_data@cons_id_map[[as.character(con_id)]] <- con_id
  }

  inverse_data@constraints <- ordered_cons
  # Batch expressions together, then split apart.
  expr_list <- list()
  for(c in ordered_cons) {
    for(arg in c@args)
      expr_list <- c(expr_list, arg)
  }
  params_to_Ab <- affine(extractor, expr_list)

  inverse_data@minimize <- inherits(problem@objective, "Minimize")
  new_prob <- ParamQuadProg(params_to_P, params_to_q, flattened_variable, params_to_Ab, variables(problem),
                            inverse_data@var_offsets, ordered_cons, parameters(problem), inverse_data@param_id_map)
  return(list(new_prob, inverse_data))
})

setMethod("invert", signature(object = "QpMatrixStuffing", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  # Retrieves the solution to the original problem.
  var_map <- inverse_data@var_offsets
  # Flip sign of opt val if maximize.
  opt_val <- solution@opt_val
  if(!(solution@status %in% ERROR_STATUS) && !inverse_data@minimize)
    opt_val <- -solution@opt_val

  primal_vars <- list()
  dual_vars <- list()
  if(!(solution@status %in% SOLUTION_PRESENT))
    return(Solution(solution@status, opt_val, primal_vars, dual_vars, solution@attr))

  # Split vectorized variable into components.
  x_opt <- values(solution@primal_vars)[[1]]
  for(var_id in names(var_map)) {
    offset <- var_map[[var_id]]
    shape <- inverse_data@var_shapes[[var_id]]
    size <- as.integer(base::prod(shape))
    primal_vars[[var_id]] <- matrix(x_opt[offset:(offset + size - 1)], nrow = shape[1], ncol = shape[2], byrow = FALSE)
  }

  if(!is.null(solution@dual_vars) && length(solution@dual_vars) > 0) {
    # Giant dual variable.
    dual_var <- values(solution@dual_vars)[[1]]
    offset <- 0
    for(constr in inverse_data@constraints) {
      # QP constraints can only have one argument.
      arg_shape <- dim(constr@args[[1]])
      dual_vars[[as.character(constr@id)]] <- matrix(dual_var[offset:(offset + size(constr@args[[1]]) - 1)], nrow = arg_shape[1], ncol = arg_shape[2], byrow = FALSE)
      offset <- offset + size(constr)
    }
  }

  return(Solution(solution@status, opt_val, primal_vars, dual_vars, solution@attr))
})


# Atom canonicalizers
Qp2QuadForm.huber_canon <- function(expr, args) {
  M <- expr@M
  x <- args[[1]]
  expr_dim <- dim(expr)
  # n <- Variable(expr_dim)
  # s <- Variable(expr_dim)
  n <- new("Variable", dim = expr_dim)
  s <- new("Variable", dim = expr_dim)

  # n^2 + 2*M*|s|
  # TODO: Make use of recursion inherent to canonicalization process and just return a power / abs expression for readability's sake.
  power_expr <- Power(n, 2)
  canon <- Qp2QuadForm.power_canon(power_expr, power_expr@args)
  n2 <- canon[[1]]
  constr_sq <- canon[[2]]

  abs_expr <- abs(s)
  canon <- EliminatePwl.abs_canon(abs_expr, abs_expr@args)
  abs_s <- canon[[1]]
  constr_abs <- canon[[2]]
  obj <- n2 + 2*M*abs_s

  constraints <- c(constr_sq, constr_abs)
  constraints <- c(constraints, list(x == s + n))
  return(list(obj, constraints))
}

Qp2QuadForm.power_canon <- function(expr, args) {
  affine_expr <- args[[1]]
  if(is.numeric(expr@p))
    p <- expr@p
  else
    p <- value(expr@p)

  if(is_constant(expr))
    return(list(Constant(value(expr)), list()))
  else if(p == 0)
    return(list(matrix(1, nrow = nrow(affine_expr), ncol = ncol(affine_expr)), list()))
  else if(p == 1)
    return(list(affine_expr, list()))
  else if(p == 2) {
    if(is(affine_expr, "Variable")) {
      affine_expr_size <- size(affine_expr)
      speye <- sparseMatrix(1:affine_expr_size, 1:affine_expr_size, x = rep(1, affine_expr_size))
      return(list(SymbolicQuadForm(affine_expr, speye, expr), list()))
    } else {
      # t <- Variable(dim(affine_expr))
      t <- new("Variable", dim = dim(affine_expr))
      t_size <- size(t)
      speye <- sparseMatrix(1:t_size, 1:t_size, x = rep(1, t_size))
      return(list(SymbolicQuadForm(t, speye, expr), list(affine_expr == t)))
    }
  }
  stop("Non-constant quadratic forms cannot be raised to a power greater than 2.")
}

Qp2QuadForm.quad_form_canon <- function(expr, args) {
  affine_expr <- expr@args[[1]]
  P <- expr@args[[2]]
  if(is(affine_expr, "Variable"))
    return(list(SymbolicQuadForm(affine_expr, P, expr), list()))
  else {
    # t <- Variable(dim(affine_expr))
    t <- new("Variable", dim = dim(affine_expr))
    return(list(SymbolicQuadForm(t, P, expr), list(affine_expr == t)))
  }
}

Qp2QuadForm.quad_over_lin_canon <- function(expr, args) {
  affine_expr <- args[[1]]
  y <- args[[2]]
  # Simplify if y has no parameters
  if(length(parameters(y)) == 0) {
    affine_expr_size <- size(affine_expr)
    speye <- sparseMatrix(1:affine_expr_size, 1:affine_expr_size, x = rep(1, affine_expr_size))
    quad_mat <- speye/value(y)
  } else {
    # TODO: This code path produces an intermediate dense matrix, but it should be sparse the whole time.
    affine_expr_size <- size(affine_expr)
    speye <- sparseMatrix(1:affine_expr_size, 1:affine_expr_size, x = rep(1, affine_expr_size))
    quad_mat <- speye/y
  }

  if(is(affine_expr, "Variable"))
    return(list(SymbolicQuadForm(affine_expr, quad_mat, expr), list()))
  else {
    # t <- Variable(dim(affine_expr))
    t <- new("Variable", dim = dim(affine_expr))
    return(list(SymbolicQuadForm(t, quad_mat, expr), list(affine_expr == t)))
  }
}

# TODO: Remove pwl canonicalize methods, use EliminatePwl reduction instead.

# Conic canonicalization methods.
Qp2QuadForm.CANON_METHODS <- list(
  Abs = Dcp2Cone.CANON_METHODS$Abs,
  CumSum = Dcp2Cone.CANON_METHODS$CumSum,
  MaxElemwise = Dcp2Cone.CANON_METHODS$MaxElemwise,
  MinElemwise = Dcp2Cone.CANON_METHODS$MinElemwise,
  SumLargest = Dcp2Cone.CANON_METHODS$SumLargest,
  MaxEntries = Dcp2Cone.CANON_METHODS$MaxEntries,
  MinEntries = Dcp2Cone.CANON_METHODS$MinEntries,
  Norm1 = Dcp2Cone.CANON_METHODS$Norm1,
  NormInf = Dcp2Cone.CANON_METHODS$NormInf,
  Indicator = Dcp2Cone.CANON_METHODS$Indicator,
  SpecialIndex = Dcp2Cone.CANON_METHODS$SpecialIndex)

# Canonicalizations that are different for QPs.
Qp2QuadForm.QUAD_CANON_METHODS <- list(
  QuadOverLin = Qp2QuadForm.quad_over_lin_canon,
  Power = Qp2QuadForm.power_canon,
  Huber = Qp2QuadForm.huber_canon,
  QuadForm = Qp2QuadForm.quad_form_canon)

Qp2QuadForm.CANON_METHODS <- c(Qp2QuadForm.CANON_METHODS, Qp2QuadForm.QUAD_CANON_METHODS)
