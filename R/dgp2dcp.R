Dgp2Dcp <- setClass("Dgp2Dcp", contains = "Canonicalization")

setMethod("accepts", signature(object = "Dgp2Dcp", problem = "Problem"), function(object, problem) {
  return(is_dgp(problem))
})

setMethod("perform", signature(object = "Dgp2Dcp", problem = "Problem"), function(object, problem) {
  if(!accepts(object, problem))
    stop("The supplied problem is not DGP")
  object@canon_methods <- DgpCanonMethods()
  tmp <- callNextMethod(object, problem)
  equiv_problem <- tmp[[1]]
  inverse_data <- tmp[[2]]
  inverse_data@problem <- problem
  return(list(object, equiv_problem, inverse_data))
})

setMethod("canonicalize_expr", signature(object = "Dgp2Dcp", expr = "Expression"), function(object, expr, args) {
  if(class(expr) %in% names(object@canon_methods))
    return(object@canon_methods[[class(expr)]](expr, args))
  else
    return(list(copy(expr, args), list()))
})

setMethod("invert", signature(object = "Dgp2Dcp", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  solution <- callNextMethod(object, solution, inverse_data)
  if(solution@status == SOLVER_ERROR)
    return(solution)
  for(vid in names(solution@primal_vars))
    solution@primal_vars[vid] <- exp(solution@primal_vars[vid])
  # f(x) = e^{F(u)}
  solution@opt_val <- exp(solution@opt_val)
  return(solution)
})

# Atom canonicalizers
# TODO: Implement sum_largest/sum_smallest.
Dgp2Dcp.add_canon <- function(expr, args) {
  if(is_scalar(expr))
    return(list(log_sum_exp(do.call("HStack", args)), list()))
  
  rows <- c()
  summands <- sapply(args, function(s) { if(is_scalar(s)) promote(s, dim(expr)) else s })
  if(length(dim(expr)) == 1) {
    for(i in 1:nrow(expr)) {
      summand_args <- lapply(summands, function(summand) { summand[i] })
      rows <- c(rows, log_sum_exp(do.call("HStack", summand_args)))
    }
    return(list(reshape_expr(bmat(rows), dim(expr)), list()))
  } else {
    for(i in 1:nrow(expr)) {
      row <- c()
      for(j in 1:ncol(expr)) {
        summand_args <- lapply(summands, function(summand) { summand[i,j] })
        rows <- c(rows, log_sum_exp(do.call("HStack", summand_args)))
      }
    }
    return(list(reshape_expr(bmat(rows), dim(expr)), list()))
  }
}

Dgp2Dcp.constant_canon <- function(expr, args) {
  # args <- list()
  return(list(Constant(log(value(expr))), list()))
}

Dgp2Dcp.div_canon <- function(expr, args) {
  # expr <- NULL
  # x / y == x * y^(-1)
  return(list(args[[1]] - args[[2]], list()))
}

Dgp2Dcp.exp_canon <- function(expr, args) {
  # expr <- NULL
  return(list(Exp(args[[1]]), list()))
}

Dgp2Dcp.eye_minus_inv_canon <- function(expr, args) {
  X <- args[[1]]
  # (I - X)^(-1) <= T iff there exists 0 <= Y <= T s.t. YX + Y <= Y.
  # Y represents log(Y) here, hence no positivity constraint.
  Y <- Variable(dim(X))
  prod <- MulExpression(Y, X)
  lhs <- Dgp2Dcp.mulexpression_canon(prod, prod@args)[[1]]
  lhs <- lhs + diag(1, nrow(prod))
  return(list(Y, list(lhs <= Y)))
}

Dgp2Dcp.geo_mean_canon <- function(expr, args) {
  out <- 0.0
  for(i in 1:length(args[[1]])) {
    x_i <- args[[1]][[i]]
    p_i <- expr@p[i]
    out <- out + p_i * x_i
  }
  return(list((1 / sum(expr@p))*out, list()))
}

Dgp2Dcp.log_canon <- function(expr, args) {
  return(list(Log(args[[1]]), list()))
}

Dgp2Dcp.mul_canon <- function(expr, args) {
  # expr <- NULL
  return(list(AddExpression(args), list()))
}

Dgp2Dcp.mulexpression_canon <- function(expr, args) {
  lhs <- args[[1]]
  rhs <- args[[2]]
  dims <- mul_dims_promote(dim(lhs), dim(rhs))
  lhs_dim <- dims[[1]]
  rhs_dim <- dims[[2]]
  lhs <- reshape_expr(lhs, lhs_dim)
  rhs <- reshape_expr(rhs, rhs_dim)
  rows <- c()
  
  # TODO: Parallelize this for large matrices.
  for(i in 1:nrow(lhs)) {
    row <- c()
    for(j in 1:ncol(rhs)) {
      hstack_args <- lapply(1:ncol(lhs), function(k) { lhs[i,k] + rhs[k,j] })
      row <- c(row, log_sum_exp(do.call("HStack", hstack_args)))
    }
    rows <- c(rows, row)
  }
  mat <- bmat(rows)
  if(!all(dim(mat) == dim(expr)))
    mat <- reshape_expr(mat, dim(expr))
  return(list(mat, list()))
}

Dgp2Dcp.nonpos_constr_canon <- function(expr, args) {
  if(length(args) != 2)
    stop("Must have exactly 2 arguments")
  return(list(NonPosConstraint(args[[1]] - args[[2]], id = id(expr)), list()))
}

Dgp2Dcp.norm1_canon <- function(expr, args) {
  if(length(args) != 1)
    stop("Must have exactly 1 argument")
  tmp <- SumEntries(args[[1]], axis = expr@axis, keepdims = expr@keepdims)
  return(Dgp2Dcp.sum_canon(tmp, tmp@args))
}

Dgp2Dcp.norm_inf_canon <- function(expr, args) {
  if(length(args) != 1)
    stop("Must have exactly 1 argument")
  tmp <- MaxEntries(args[[1]], axis = expr@axis, keepdims = expr@keepdims)
  return(EliminatePwl.max_entries_canon(tmp, tmp@args))
}

Dgp2Dcp.one_minus_pos_canon <- function(expr, args) {
  return(list(Log(expr@ones - Exp(args[[1]])), list()))
}

Dgp2Dcp.parameter_canon <- function(expr, args) {
  # args <- list()
  return(list(Parameter(log(value(expr)), name = name(expr)), list()))
}

Dgp2Dcp.pf_eigenvalue_canon <- function(expr, args) {
  X <- args[[1]]
  # rho(X) <= lambda iff there exists v s.t. Xv <= lambda v.
  # v and lambda represent log variables, hence no positivity constraints.
  lambd <- Variable()
  v <- Variable(nrow(X))
  lhs <- MulExpression(X, v)
  rhs <- lambd*v
  lhs <- Dgp2Dcp.mulexpression_canon(lhs, lhs@args)[[1]]
  rhs <- Dgp2Dcp.mul_canon(rhs, rhs@args)[[1]]
  return(list(lambd, list(lhs <= rhs)))
}

Dgp2Dcp.pnorm_canon <- function(expr, args) {
  x <- args[[1]]
  p <- expr@original_p
  if(is.null(dim(x)))
    x <- promote(p, c(1))
  if(is.na(expr@axis) || length(dim(x)) == 1) {
    x <- Vec(x)
    hstack_args <- lapply(x, function(xi) { xi^p })
    return(list((1.0/p) * log_sum_exp(do.call("HStack", hstack_args)), list()))
  }
  
  if(expr@axis == 2)
    x <- t(x)
  
  rows <- c()
  for(i in 1:nrow(x)) {
    row <- x[i]
    hstack_args <- lapply(row, function(xi) { xi^p })
    rows <- c(rows, (1.0/p)*log_sum_exp(do.call("HStack", hstack_args)))
  }
  return(list(vstack(rows), list()))
}

Dgp2Dcp.power_canon <- function(expr, args) {
  # y = log(x); x^p --> exp(y^p) --> p*log(exp(y)) = p*y.
  return(list(expr@p*args[[1]], list()))
}

Dgp2Dcp.prod_canon <- function(expr, args) {
  return(list(SumEntries(args[[1]], axis = expr@axis, keepdims = expr@keepdims), list()))
}

Dgp2Dcp.quad_form_canon <- function(expr, args) {
  x <- args[[1]]
  P <- args[[2]]
  elems <- list()
  for(i in 1:nrow(P)) {
    for(j in 1:nrow(P))
      elems <- c(elems, P[i,j] + x[i] + x[j])
  }
  return(list(log_sum_exp(do.call("HStack", elems)), list()))
}

Dgp2Dcp.quad_over_lin_canon <- function(expr, args) {
  x <- Vec(args[[1]])
  y <- args[[2]]
  numerator <- sum(sapply(x, function(xi) { 2*xi }))
  return(list(numerator - y, list()))
}

Dgp2Dcp.sum_canon <- function(expr, args) {
  X <- args[[1]]
  if(is.na(expr@axis)) {
    x <- Vec(X)
    summation <- do.call("sum", args = lapply(x, function(xi) { xi }))
    canon <- Dgp2Dcp.add_canon(summation, summation@args)[[1]]
    return(list(reshape_expr(canon, dim(expr)), list()))
  }
  
  if(expr@axis == 2)
    X <- t(X)
  
  rows <- list()
  for(i in 1:nrow(X)) {
    x <- Vec(X[i])
    summation <- do.call("sum", args = lapply(x, function(xi) { xi }))
    canon <- Dgp2Dcp.add_canon(summation, summation@args)[[1]]
    rows <- c(rows, canon)
  }
  canon <- do.call("HStack", rows)
  return(list(reshape_expr(canon, dim(expr)), list()))
}

Dgp2Dcp.trace_canon <- function(expr, args) {
  diag_sum <- sum(Diag(args[[1]]))
  return(Dgp2Dcp.add_canon(diag_sum, diag_sum@args))
}

Dgp2Dcp.zero_constr_canon <- function(expr, args) {
  if(length(args) != 2)
    stop("Must have exactly 2 arguments")
  return(list(ZeroConstraint(args[[1]] - args[[2]], id = id(expr)), list()))
}

Dgp2Dcp.CANON_METHODS <- list(AddExpression = Dgp2Dcp.add_canon,
                              Constant = Dgp2Dcp.constant_canon,
                              DivExpression = Dgp2Dcp.div_canon,
                              Exp = Dgp2Dcp.exp_canon,
                              EyeMinusInv = Dgp2Dcp.eye_minus_inv_canon,
                              GeoMean = Dgp2Dcp.geo_mean_canon,
                              Log = Dgp2Dcp.log_canon,
                              MulExpression = Dgp2Dcp.mulexpression_canon,
                              Multiply = Dgp2Dcp.mul_canon,
                              Norm1 = Dgp2Dcp.norm1_canon,
                              NormInf = Dgp2Dcp.norm_inf_canon,
                              OneMinusPos = Dgp2Dcp.one_minus_pos_canon,
                              Parameter = Dgp2Dcp.parameter_canon,
                              PfEigenvalue = Dgp2Dcp.pf_eigenvalue_canon,
                              Pnorm = Dgp2Dcp.pnorm_canon,
                              Power = Dgp2Dcp.power_canon,
                              Prod = Dgp2Dcp.prod_canon,
                              QuadForm = Dgp2Dcp.quad_form_canon,
                              QuadOverLin = Dgp2Dcp.quad_over_lin_canon,
                              Trace = Dgp2Dcp.trace_canon,
                              SumEntries = Dgp2Dcp.sum_canon,
                              Variable = NULL,
                              
                              MaxEntries = EliminatePwl.CANON_METHODS$MaxEntries,
                              MinEntries = EliminatePwl.CANON_METHODS$MinEntries,
                              MaxElemwise = EliminatePwl.CANON_METHODS$MaxElemwise,
                              MinElemwise = EliminatePwl.CANON_METHODS$MinElemwise)

# Canonicalization of DGPs is a stateful procedure, hence the need for a class.
.DgpCanonMethods <- setClass("DgpCanonMethods", representation(.variables = "list"), prototype(.variables = list()), contains = "list")
DgpCanonMethods <- function(...) { .DgpCanonMethods(...) }

setMethod("[", signature(x = "DgpCanonMethods", i = "character", j = "missing", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  if(i == "Variable") {
    # TODO: Check scoping of x here is correct.
    variable_canon <- function(variable, args) {
      # Swap out positive variables for unconstrained variables.
      if(id(variable) %in% names(x@.variables))
        return(list(variable, list()))
      else {
        log_variable <- Variable(dim(variable), var_id = id(variable))
        x@.variables[as.character(id(variable))] <- log_variable
        return(list(log_variable), list())
      }
    }
    return(variable_canon)
  } else
    return(Dgp2Dcp.CANON_METHODS[[i]])
})
