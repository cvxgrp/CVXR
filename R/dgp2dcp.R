#'
#' Reduce DGP problems to DCP problems.
#' 
#' This reduction takes as input a DGP problem and returns an equivalent DCP
#' problem. Because every (generalized) geometric program is a DGP problem,
#' this reduction can be used to convert geometric programs into convex form.
#' @rdname Dgp2Dcp-class
Dgp2Dcp <- setClass("Dgp2Dcp", contains = "Canonicalization")

#' @describeIn Dgp2Dcp Is the problem DGP?
setMethod("accepts", signature(object = "Dgp2Dcp", problem = "Problem"), function(object, problem) {
  return(is_dgp(problem))
})

#' @param object A \linkS4class{Dgp2Dcp} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn Dgp2Dcp Converts the DGP problem to a DCP problem.
setMethod("perform", signature(object = "Dgp2Dcp", problem = "Problem"), function(object, problem) {
  if(!accepts(object, problem))
    stop("The supplied problem is not DGP")
  object@canon_methods <- DgpCanonMethods()
  tmp <- callNextMethod(object, problem)
  object <- tmp[[1]]
  equiv_problem <- tmp[[2]]
  inverse_data <- tmp[[3]]
  inverse_data@problem <- problem
  return(list(object, equiv_problem, inverse_data))
})

#' @param expr An \linkS4class{Expression} object corresponding to the DGP problem.
#' @param args A list of values corresponding to the DGP expression
#' @describeIn Dgp2Dcp Canonicalizes each atom within an Dgp2Dcp expression.
setMethod("canonicalize_expr", "Dgp2Dcp", function(object, expr, args) {
  if(class(expr) %in% names(object@canon_methods))
    return(object@canon_methods[[class(expr)]](expr, args))
  else
    return(list(copy(expr, args), list()))
})

#' @param solution A \linkS4class{Solution} object to invert.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn Dgp2Dcp Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "Dgp2Dcp", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  solution <- callNextMethod(object, solution, inverse_data)
  if(solution@status == SOLVER_ERROR)
    return(solution)
  for(vid in names(solution@primal_vars))
    solution@primal_vars[[vid]] <- exp(solution@primal_vars[[vid]])
  # f(x) = e^{F(u)}
  solution@opt_val <- exp(solution@opt_val)
  return(solution)
})

# Atom canonicalizers
# TODO: Implement sum_largest/sum_smallest.
#' 
#' Dgp2Dcp canonicalizer for the addition atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the addition atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.add_canon <- function(expr, args) {
  if(is_scalar(expr))
    return(list(log_sum_exp(do.call("HStack", args)), list()))
  
  rows <- list()
  summands <- lapply(args, function(s) { if(is_scalar(s)) promote(s, dim(expr)) else s })
  if(length(dim(expr)) == 1) {
    for(i in 1:nrow(expr)) {
      summand_args <- lapply(summands, function(summand) { summand[i] })
      row <- log_sum_exp(do.call("HStack", summand_args))
      rows <- c(rows, list(row))
    }
    return(list(reshape_expr(bmat(rows), dim(expr)), list()))
  } else {
    for(i in 1:nrow(expr)) {
      row <- list()
      for(j in 1:ncol(expr)) {
        summand_args <- lapply(summands, function(summand) { summand[i,j] })
        row <- c(row, list(log_sum_exp(do.call("HStack", summand_args))))
      }
      rows <- c(rows, list(row))
    }
    return(list(reshape_expr(bmat(rows), dim(expr)), list()))
  }
}

#' 
#' Dgp2Dcp canonicalizer for the constant atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the constant atom of a DGP expression, 
#' where the returned expression is the DCP equivalent resulting 
#' from the log of the expression.
Dgp2Dcp.constant_canon <- function(expr, args) {
  # args <- list()
  return(list(Constant(log(value(expr))), list()))
}

#' 
#' Dgp2Dcp canonicalizer for the division atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the division atom of a DGP expression, 
#' where the returned expression is the log transformed DCP equivalent.
Dgp2Dcp.div_canon <- function(expr, args) {
  # expr <- NULL
  # x / y == x * y^(-1)
  return(list(args[[1]] - args[[2]], list()))
}

#' 
#' Dgp2Dcp canonicalizer for the exp atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the exp atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.exp_canon <- function(expr, args) {
  # expr <- NULL
  return(list(Exp(args[[1]]), list()))
}

#' 
#' Dgp2Dcp canonicalizer for the (I - X)^{-1} atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the (I - X)^{-1} atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.eye_minus_inv_canon <- function(expr, args) {
  X <- args[[1]]
  # (I - X)^(-1) <= T iff there exists 0 <= Y <= T s.t. YX + Y <= Y.
  # Y represents log(Y) here, hence no positivity constraint.
  # Y <- Variable(dim(X))
  Y <- new("Variable", dim = dim(X))
  prod <- MulExpression(Y, X)
  lhs <- Dgp2Dcp.mulexpression_canon(prod, prod@args)[[1]]
  lhs <- lhs + diag(1, nrow(prod))
  return(list(Y, list(lhs <= Y)))
}

#' 
#' Dgp2Dcp canonicalizer for the geometric mean atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the geometric mean atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.geo_mean_canon <- function(expr, args) {
  out <- 0.0
  for(i in seq_along(args[[1]])) {
    x_i <- args[[1]][i]
    p_i <- expr@p[i]
    out <- out + p_i * x_i
  }
  return(list((1 / sum(expr@p))*out, list()))
}

#' 
#' Dgp2Dcp canonicalizer for the log atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the log atom of a DGP expression,
#' where the returned expression is the log of the original expression..
Dgp2Dcp.log_canon <- function(expr, args) {
  return(list(Log(args[[1]]), list()))
}

#' 
#' Dgp2Dcp canonicalizer for the multiplication atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the multiplication atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.mul_canon <- function(expr, args) {
  # expr <- NULL
  return(list(AddExpression(args), list()))
}

#' 
#' Dgp2Dcp canonicalizer for the multiplication expression atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the multiplication expression atom 
#' of a DGP expression, where the returned expression is the transformed 
#' DCP equivalent.
Dgp2Dcp.mulexpression_canon <- function(expr, args) {
  lhs <- args[[1]]
  rhs <- args[[2]]
  dims <- mul_dims_promote(dim(lhs), dim(rhs))
  lhs_dim <- dims[[1]]
  rhs_dim <- dims[[2]]
  lhs <- reshape_expr(lhs, lhs_dim)
  rhs <- reshape_expr(rhs, rhs_dim)
  rows <- list()
  
  # TODO: Parallelize this for large matrices.
  for(i in 1:nrow(lhs)) {
    row <- list()
    for(j in 1:ncol(rhs)) {
      hstack_args <- lapply(1:ncol(lhs), function(k) { lhs[i,k] + rhs[k,j] })
      row <- c(row, list(log_sum_exp(do.call("HStack", hstack_args))))
    }
    rows <- c(rows, list(row))
  }
  mat <- bmat(rows)
  if(!all(dim(mat) == dim(expr)))
    mat <- reshape_expr(mat, dim(expr))
  return(list(mat, list()))
}

#' 
#' Dgp2Dcp canonicalizer for the non-positive constraint atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the non-positive contraint atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.nonpos_constr_canon <- function(expr, args) {
  if(length(args) != 2)
    stop("Must have exactly 2 arguments")
  return(list(NonPosConstraint(args[[1]] - args[[2]], id = id(expr)), list()))
}

#' 
#' Dgp2Dcp canonicalizer for the 1 norm atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the norm1 atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.norm1_canon <- function(expr, args) {
  if(length(args) != 1)
    stop("Must have exactly 1 argument")
  tmp <- SumEntries(args[[1]], axis = expr@axis, keepdims = expr@keepdims)
  return(Dgp2Dcp.sum_canon(tmp, tmp@args))
}

#' 
#' Dgp2Dcp canonicalizer for the infinite norm atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the infinity norm atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.norm_inf_canon <- function(expr, args) {
  if(length(args) != 1)
    stop("Must have exactly 1 argument")
  tmp <- MaxEntries(args[[1]], axis = expr@axis, keepdims = expr@keepdims)
  return(EliminatePwl.max_entries_canon(tmp, tmp@args))
}

#' 
#' Dgp2Dcp canonicalizer for the 1-x atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the 1-x with 0 < x < 1 atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.one_minus_pos_canon <- function(expr, args) {
  return(list(Log(expr@.ones - Exp(args[[1]])), list()))
}

#' 
#' Dgp2Dcp canonicalizer for the parameter atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the parameter atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.parameter_canon <- function(expr, args) {
  # args <- list()
  return(list(Parameter(log(value(expr)), name = name(expr)), list()))
}

#' 
#' Dgp2Dcp canonicalizer for the spectral radius atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the spectral radius atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
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

#' 
#' Dgp2Dcp canonicalizer for the p norm atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the pnorm atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.pnorm_canon <- function(expr, args) {
  x <- args[[1]]
  p <- expr@original_p
  if(is.null(dim(x)))
    x <- promote(p, c(1))
  if(is.na(expr@axis) || length(dim(x)) == 1) {
    x <- Vec(x)
    # hstack_args <- lapply(seq_len(size(x)), function(j) { x[j]^p })
    # return(list((1.0/p) * log_sum_exp(do.call("HStack", hstack_args)), list()))
    return(list((1.0/p) * log_sum_exp(x^p), list()))
  }
  
  if(expr@axis == 2)
    x <- t(x)
  
  rows <- list()
  for(i in 1:nrow(x)) {
    row <- x[i,]
    # hstack_args <- lapply(seq_len(size(row)), function(j) { row[j]^p })
    # rows <- c(rows, list((1.0/p)*log_sum_exp(do.call("HStack", hstack_args))))
    rows <- c(rows, list((1.0/p) * log_sum_exp(row^p)))
  }
  return(list(do.call("VStack", rows), list()))
}

#' 
#' Dgp2Dcp canonicalizer for the power atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the power atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.power_canon <- function(expr, args) {
  # y = log(x); x^p --> exp(y^p) --> p*log(exp(y)) = p*y.
  return(list(expr@p*args[[1]], list()))
}

#' 
#' Dgp2Dcp canonicalizer for the product atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the product atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.prod_canon <- function(expr, args) {
  return(list(SumEntries(args[[1]], axis = expr@axis, keepdims = expr@keepdims), list()))
}

#' 
#' Dgp2Dcp canonicalizer for the quadratic form atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the quadratic form atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.quad_form_canon <- function(expr, args) {
  x <- args[[1]]
  P <- args[[2]]
  elems <- list()
  for(i in 1:nrow(P)) {
    for(j in 1:nrow(P))
      elems <- c(elems, list(P[i,j] + x[i] + x[j]))
  }
  return(list(log_sum_exp(do.call("HStack", elems)), list()))
}

#' 
#' Dgp2Dcp canonicalizer for the quadratic over linear term atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the quadratic over linear atom of a 
#' DGP expression, where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.quad_over_lin_canon <- function(expr, args) {
  x <- Vec(args[[1]])
  y <- args[[2]]
  numerator <- 2*sum(x)
  return(list(numerator - y, list()))
}

#' 
#' Dgp2Dcp canonicalizer for the sum atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the sum atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.sum_canon <- function(expr, args) {
  X <- args[[1]]
  new_dim <- dim(expr)
  if(is.null(new_dim))
    new_dim <- c(1,1)
  
  if(is.na(expr@axis)) {
    x <- Vec(X)
    summation <- sum(x)
    canon <- Dgp2Dcp.add_canon(summation, summation@args)[[1]]
    return(list(reshape_expr(canon, new_dim), list()))
  }
  
  if(expr@axis == 2)
    X <- t(X)
  
  rows <- list()
  for(i in 1:nrow(X)) {
    x <- Vec(X[i,])
    summation <- sum(x)
    canon <- Dgp2Dcp.add_canon(summation, summation@args)[[1]]
    rows <- c(rows, list(canon))
  }
  canon <- do.call("HStack", rows)
  return(list(reshape_expr(canon, new_dim), list()))
}

#' 
#' Dgp2Dcp canonicalizer for the trace atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the trace atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.trace_canon <- function(expr, args) {
  diag_sum <- sum(Diag(args[[1]]))
  return(Dgp2Dcp.add_canon(diag_sum, diag_sum@args))
}

#' 
#' Dgp2Dcp canonicalizer for the zero constraint atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the zero constraint atom of a DGP expression, 
#' where the returned expression is the transformed DCP equivalent.
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
                              ProdEntries = Dgp2Dcp.prod_canon,
                              QuadForm = Dgp2Dcp.quad_form_canon,
                              QuadOverLin = Dgp2Dcp.quad_over_lin_canon,
                              Trace = Dgp2Dcp.trace_canon,
                              SumEntries = Dgp2Dcp.sum_canon,
                              Variable = NULL,
                              
                              MaxEntries = EliminatePwl.CANON_METHODS$MaxEntries,
                              MinEntries = EliminatePwl.CANON_METHODS$MinEntries,
                              MaxElemwise = EliminatePwl.CANON_METHODS$MaxElemwise,
                              MinElemwise = EliminatePwl.CANON_METHODS$MinElemwise)

#' 
#' DGP canonical methods class.
#' 
#' Canonicalization of DGPs is a stateful procedure, hence the need for a class.
#' 
#' @rdname DgpCanonMethods-class
.DgpCanonMethods <- setClass("DgpCanonMethods", representation(.variables = "list"), prototype(.variables = list()), contains = "list")
DgpCanonMethods <- function(...) { .DgpCanonMethods(...) }

#' @param x A \linkS4class{DgpCanonMethods} object.
#' @describeIn DgpCanonMethods Returns the name of all the canonicalization methods
setMethod("names", signature(x = "DgpCanonMethods"), function(x) { names(Dgp2Dcp.CANON_METHODS) })

# TODO: How to implement this with S4 setMethod? Signature is x = "DgpCanonMethods", i = "character", j = "missing".
'[[.DgpCanonMethods' <- function(x, i, j, ..., exact = TRUE) { do.call("$", list(x, i)) }

#' @param name The name of the atom or expression to canonicalize.
#' @describeIn DgpCanonMethods Returns either a canonicalized variable or 
#'  a corresponding Dgp2Dcp canonicalization method 
setMethod("$", signature(x = "DgpCanonMethods"), function(x, name) {
  if(name == "Variable") {
    # TODO: Check scoping of x here is correct.
    variable_canon <- function(variable, args) {
      args <- NULL
      vid <- as.character(variable@id)
      # Swap out positive variables for unconstrained variables.
      if(vid %in% names(x@.variables))
        return(list(x@.variables[[vid]], list()))
      else {
        # log_variable <- Variable(dim(variable), id = variable@id)
        log_variable <- new("Variable", dim = dim(variable), id = variable@id)
        x@.variables[[vid]] <- log_variable
        return(list(log_variable, list()))
      }
    }
    return(variable_canon)
  } else
    return(Dgp2Dcp.CANON_METHODS[[name]])
})
