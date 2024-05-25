#'
#' The InverseData class.
#'
#' This class represents the data encoding an optimization problem.
#'
#' @rdname InverseData-class
.InverseData <- setClass("InverseData", representation(problem = "Problem", id_map = "list", var_offsets = "list", x_length = "numeric", var_dims = "list",
                                                       id2var = "list", real2imag = "list", id2cons = "list", cons_id_map = "list", r = "numeric", minimize = "logical", 
                                                       sorted_constraints = "list", is_mip = "logical"),
                         prototype(id_map = list(), var_offsets = list(), x_length = NA_real_, var_dims = list(), id2var = list(),
                                   real2imag = list(), id2cons = list(), cons_id_map = list(), r = NA_real_, minimize = NA, sorted_constraints = list(), is_mip = NA))

InverseData <- function(problem) { .InverseData(problem = problem) }

setMethod("initialize", "InverseData", function(.Object, ..., problem, id_map = list(), var_offsets = list(), x_length = NA_real_, var_dims = list(), id2var = list(), real2imag = list(), id2cons = list(), cons_id_map = list(), r = NA_real_, minimize = NA, sorted_constraints = list(), is_mip = NA) {
  # Basic variable offset information
  varis <- variables(problem)
  varoffs <- get_var_offsets(.Object, varis)
  .Object@id_map <- varoffs$id_map
  .Object@var_offsets <- varoffs$var_offsets
  .Object@x_length <- varoffs$x_length
  .Object@var_dims <- varoffs$var_dims
  
  # Map of variable id to variable
  .Object@id2var <- stats::setNames(varis, sapply(varis, function(var) { as.character(id(var)) }))
  
  # Map of real to imaginary parts of complex variables
  var_comp <- Filter(is_complex, varis)
  .Object@real2imag <- lapply(seq_along(var_comp), function(i) { get_id() })
  names(.Object@real2imag) <- sapply(var_comp, function(var) { as.character(id(var)) })
  constrs <- constraints(problem)
  constr_comp <- Filter(is_complex, constrs)
  constr_dict <- lapply(seq_along(constr_comp), function(i) { get_id() })
  names(constr_dict) <- sapply(constr_comp, function(cons) { as.character(id(cons)) })
  .Object@real2imag <- utils::modifyList(.Object@real2imag, constr_dict)
  
  # Map of constraint id to constraint
  .Object@id2cons <- stats::setNames(constrs, sapply(constrs, function(cons) { as.character(id(cons)) }))
  .Object@cons_id_map <- list()
  .Object@r <- r
  .Object@minimize <- minimize
  .Object@sorted_constraints <- sorted_constraints
  .Object@is_mip <- is_mip
  return(.Object)
})

setMethod("get_var_offsets", signature(object = "InverseData", variables = "list"), function(object, variables) {
  var_dims <- list()
  var_offsets <- list()
  id_map <- list()
  vert_offset <- 0
  for(x in variables) {
    var_dims[[as.character(id(x))]] <- dim(x)
    var_offsets[[as.character(id(x))]] <- vert_offset
    id_map[[as.character(id(x))]] <- list(vert_offset, size(x))
    vert_offset <- vert_offset + size(x)
  }
  return(list(id_map = id_map, var_offsets = var_offsets, x_length = vert_offset, var_dims = var_dims))
})

# TODO: Find best format for sparse matrices.
.CoeffExtractor <- setClass("CoeffExtractor", representation(inverse_data = "InverseData", id_map = "list", N = "numeric", var_dims = "list"),
                            prototype(id_map = list(), N = NA_real_, var_dims = list()))

CoeffExtractor <- function(inverse_data) { .CoeffExtractor(inverse_data = inverse_data) }

setMethod("initialize", "CoeffExtractor", function(.Object, inverse_data, id_map, N, var_dims) {
  .Object@id_map <- inverse_data@var_offsets
  .Object@N <- inverse_data@x_length
  .Object@var_dims <- inverse_data@var_dims
  return(.Object)
})

setMethod("get_coeffs", signature(object = "CoeffExtractor", expr = "Expression"), function(object, expr) {
  if(is_constant(expr))
    return(constant(object, expr))
  else if(is_affine(expr))
    return(affine(object, expr))
  else if(is_quadratic(expr))
    return(quad_form(object, expr))
  else
    stop("Unknown expression type ", class(expr))
})

setMethod("constant", signature(object = "CoeffExtractor", expr = "Expression"), function(object, expr) {
  size <- size(expr)
  list(Matrix(nrow = size, ncol = object@N), as.vector(value(expr)))
})

setMethod("affine", signature(object = "CoeffExtractor", expr = "list"), function(object, expr) {
  if(length(expr) == 0)
    size <- 0
  else
    size <- sum(sapply(expr, size))
  op_list <- lapply(expr, function(e) { canonical_form(e)[[1]] })
  VIJb <- get_problem_matrix(op_list, object@id_map)
  V <- VIJb[[1]]
  I <- VIJb[[2]] + 1   # TODO: Convert 0-indexing to 1-indexing in get_problem_matrix.
  J <- VIJb[[3]] + 1
  b <- VIJb[[4]]
  A <- sparseMatrix(i = I, j = J, x = V, dims = c(size, object@N))
  list(A, as.vector(b))
})

setMethod("affine", signature(object = "CoeffExtractor", expr = "Expression"), function(object, expr) {
  affine(object, list(expr))
})

setMethod("extract_quadratic_coeffs", "CoeffExtractor", function(object, affine_expr, quad_forms) {
  # Assumes quadratic forms all have variable arguments. Affine expressions can be anything.
  # Extract affine data.
  affine_problem <- Problem(Minimize(affine_expr), list())
  affine_inverse_data <- InverseData(affine_problem)
  affine_id_map <- affine_inverse_data@id_map
  affine_var_dims <- affine_inverse_data@var_dims
  extractor <- CoeffExtractor(affine_inverse_data)
  cb <- affine(extractor, expr(affine_problem@objective))
  c <- cb[[1]]
  b <- cb[[2]]

  # Combine affine data with quadratic forms.
  coeffs <- list()
  for(var in variables(affine_problem)) {
    var_id <- as.character(id(var))
    if(var_id %in% names(quad_forms)) {
      orig_id <- as.character(id(quad_forms[[var_id]][[3]]@args[[1]]))
      var_offset <- affine_id_map[[var_id]][[1]]
      var_size <- affine_id_map[[var_id]][[2]]
      if(!any(is.na(value(quad_forms[[var_id]][[3]]@P)))) {
        P <- value(quad_forms[[var_id]][[3]]@P)
        if(is(P, "sparseMatrix"))
          P <- as.matrix(P)
        c_part <- c[1, (var_offset + 1):(var_offset + var_size)]
        ## P <- sweep(P, MARGIN = 2, FUN = "*", c_part)
        .Call('_CVXR_sweep_in_place', PACKAGE = 'CVXR', P, c_part)
      } else
        P <- sparseMatrix(i = 1:var_size, j = 1:var_size, x = c[1, (var_offset + 1):(var_offset + var_size)])
      if(orig_id %in% names(coeffs)) {
        coeffs[[orig_id]]$P <- coeffs[[orig_id]]$P + P
        coeffs[[orig_id]]$q <- coeffs[[orig_id]]$q + rep(0, nrow(P))
      } else {
        coeffs[[orig_id]] <- list()
        coeffs[[orig_id]]$P <- P
        coeffs[[orig_id]]$q <- rep(0, nrow(P))
      }
    } else {
      var_offset <- affine_id_map[[var_id]][[1]]
      var_size <- as.integer(prod(affine_var_dims[[var_id]]))
      if(var_id %in% names(coeffs)) {
        coeffs[[var_id]]$P <- coeffs[[var_id]]$P + sparseMatrix(i = c(), j = c(), dims = c(var_size, var_size))
        coeffs[[var_id]]$q <- coeffs[[var_id]]$q + c[1, (var_offset + 1):(var_offset + var_size)]
      } else {
        coeffs[[var_id]] <- list()
        coeffs[[var_id]]$P <- sparseMatrix(i = c(), j = c(), dims = c(var_size, var_size))
        coeffs[[var_id]]$q <- c[1, (var_offset + 1):(var_offset + var_size)]
      }
    }
  }
  return(list(coeffs, b))
})

setMethod("coeff_quad_form", signature(object = "CoeffExtractor", expr = "Expression"), function(object, expr) {
  # Extract quadratic, linear constant parts of a quadratic objective.
  # Insert no-op such that root is never a quadratic form, for easier processing.
  root <- LinOp(NO_OP, dim(expr), list(expr), list())

  # Replace quadratic forms with dummy variables.
  tmp <- replace_quad_forms(root, list())
  root <- tmp[[1]]
  quad_forms <- tmp[[2]]

  # Calculate affine parts and combine them with quadratic forms to get the coefficients.
  tmp <- extract_quadratic_coeffs(object, root$args[[1]], quad_forms)
  coeffs <- tmp[[1]]
  constant <- tmp[[2]]

  # Restore expression.
  root$args[[1]] <- restore_quad_forms(root$args[[1]], quad_forms)

  # Sort variables corresponding to their starting indices in ascending order.
  shuffle <- order(unlist(object@id_map), decreasing = FALSE)
  offsets <- object@id_map[shuffle]

  # Concatenate quadratic matrices and vectors.
  P <- Matrix(nrow = 0, ncol = 0)
  q <- c()
  for(var_id in names(offsets)) {
    offset <- offsets[[var_id]]
    if(var_id %in% names(coeffs)) {
      P <- bdiag(P, coeffs[[var_id]]$P)
      q <- c(q, coeffs[[var_id]]$q)
    } else {
      var_dim <- object@var_dims[[var_id]]
      size <- as.integer(prod(var_dim))
      P <- bdiag(P, Matrix(0, nrow = size, ncol = size))
      q <- c(q, rep(0, size))
    }
  }

  if(!(nrow(P) == ncol(P) && ncol(P) == object@N) || length(q) != object@N)
    stop("Resulting quadratic form does not have appropriate dimensions.")

  if(length(constant) != 1)
    stop("Constant must be a scalar.")
  return(list(P, q, constant[1]))
})

# Helper function for getting args of an Expression or LinOp.
get_args <- function(x) {
  if(is(x, "Expression"))
    return(x@args)
  else
    return(x$args)
}

# Helper function for setting args of an Expression or LinOp.
set_args <- function(x, idx, val) {
  if(is(x, "Expression"))
    x@args[[idx]] <- val
  else
    x$args[[idx]] <- val
  return(x)
}

replace_quad_forms <- function(expr, quad_forms) {
  nargs <- length(get_args(expr))
  if(nargs == 0)
    return(list(expr, quad_forms))
  
  for(idx in 1:nargs) {
    arg <- get_args(expr)[[idx]]
    if(is(arg, "SymbolicQuadForm") || is(arg, "QuadForm")) {
      tmp <- replace_quad_form(expr, idx, quad_forms)
      expr <- tmp[[1]]
      quad_forms <- tmp[[2]]
    } else {
      tmp <- replace_quad_forms(arg, quad_forms)
      expr <- set_args(expr, idx, tmp[[1]])
      quad_forms <- tmp[[2]]
    }
  }
  return(list(expr, quad_forms))
}

replace_quad_form <- function(expr, idx, quad_forms) {
  quad_form <- get_args(expr)[[idx]]
  # placeholder <- Variable(dim(quad_form))
  placeholder <- new("Variable", dim = dim(quad_form))
  expr <- set_args(expr, idx, placeholder)
  quad_forms[[as.character(id(placeholder))]] <- list(expr, idx, quad_form)
  return(list(expr, quad_forms))
}

restore_quad_forms <- function(expr, quad_forms) {
  nargs <- length(get_args(expr))
  if(nargs == 0)
    return(expr)
  
  for(idx in 1:nargs) {
    arg <- get_args(expr)[[idx]]
    if(is(arg, "Variable") && as.character(id(arg)) %in% names(quad_forms))
      expr <- set_args(expr, idx, quad_forms[[as.character(id(arg))]][[3]])
    else {
      arg <- restore_quad_forms(arg, quad_forms)
      expr <- set_args(expr, idx, arg)
    }
  }
  return(expr)
}

# #
# # The QuadCoeffExtractor class
# #
# # Given a quadratic expression of size m*n, this class extracts the coefficients
# # (Ps, Q, R) such that the (i,j) entry of the expression is given by
# # t(X) %*% Ps[[k]] %*% x + Q[k,] %*% x + R[k]
# # where k = i + j*m. x is the vectorized variables indexed by id_map
# #
# setClass("QuadCoeffExtractor", representation(id_map = "list", N = "numeric"))
#
# get_coeffs.QuadCoeffExtractor <- function(object, expr) {
#   if(is_constant(expr))
#     return(.coeffs_constant(object, expr))
#   else if(is_affine(expr))
#     return(.coeffs_affine(object, expr))
#   else if(is(expr, "AffineProd"))
#     return(.coeffs_affine_prod(object, expr))
#   else if(is(expr, "QuadOverLin"))
#     return(.coeffs_quad_over_lin(object, expr))
#   else if(is(expr, "Power"))
#     return(.coeffs_power(object, expr))
#   else if(is(expr, "MatrixFrac"))
#     return(.coeffs_matrix_frac(object, expr))
#   else if(is(expr, "AffAtom"))
#     return(.coeffs_affine_atom(object, expr))
#   else
#     stop("Unknown expression type: ", class(expr))
# }
#
# .coeffs_constant.QuadCoeffExtractor <- function(object, expr) {
#   if(is_scalar(expr)) {
#     sz <- 1
#     R <- matrix(value(expr))
#   } else {
#     sz <- prod(size(expr))
#     R <- matrix(value(expr), nrow = sz)    # TODO: Check if this should be transposed
#   }
#   Ps <- lapply(1:sz, function(i) { sparseMatrix(i = c(), j = c(), dims = c(object@N, object@N)) })
#   Q <- sparseMatrix(i = c(), j = c(), dims = c(sz, object@N))
#   list(Ps, Q, R)
# }
#
# .coeffs_affine.QuadCoeffExtractor <- function(object, expr) {
#   sz <- prod(size(expr))
#   canon <- canonical_form(expr)
#   prob_mat <- get_problem_matrix(list(create_eq(canon[[1]])), object@id_map)
#
#   V <- prob_mat[[1]]
#   I <- prob_mat[[2]]
#   J <- prob_mat[[3]]
#   R <- prob_mat[[4]]
#
#   Q <- sparseMatrix(i = I, j = J, x = V, dims = c(sz, object@N))
#   Ps <- lapply(1:sz, function(i) { sparseMatrix(i = c(), j = c(), dims = c(object@N, object@N)) })
#   list(Ps, Q, as.vector(R))   # TODO: Check if R is flattened correctly
# }
#
# .coeffs_affine_prod.QuadCoeffExtractor <- function(object, expr) {
#   XPQR <- .coeffs_affine(expr@args[[1]])
#   YPQR <- .coeffs_affine(expr@args[[2]])
#   Xsize <- size(expr@args[[1]])
#   Ysize <- size(expr@args[[2]])
#
#   XQ <- XPQR[[2]]
#   XR <- XPQR[[3]]
#   YQ <- YPQR[[2]]
#   YR <- YPQR[[3]]
#   m <- Xsize[1]
#   p <- Xsize[2]
#   n <- Ysize[2]
#
#   Ps  <- list()
#   Q <- sparseMatrix(i = c(), j = c(), dims = c(m*n, object@N))
#   R <- rep(0, m*n)
#
#   ind <- 0
#   for(j in 1:n) {
#     for(i in 1:m) {
#       M <- sparseMatrix(i = c(), j = c(), dims = c(object@N, object@N))
#       for(k in 1:p) {
#         Xind <- k*m + i
#         Yind <- j*p + k
#
#         a <- XQ[Xind,]
#         b <- XR[Xind]
#         c <- YQ[Yind,]
#         d <- YR[Yind]
#
#         M <- M + t(a) %*% c
#         Q[ind,] <- Q[ind,] + b*c + d*a
#         R[ind] <- R[ind] + b*d
#       }
#       Ps <- c(Ps, Matrix(M, sparse = TRUE))
#       ind <- ind + 1
#     }
#   }
#   list(Ps, Matrix(Q, sparse = TRUE), R)
# }
#
# .coeffs_quad_over_lin.QuadCoeffExtractor <- function(object, expr) {
#   coeffs <- .coeffs_affine(object, expr@args[[1]])
#   A <- coeffs[[2]]
#   b <- coeffs[[3]]
#
#   P <- t(A) %*% A
#   q <- Matrix(2*t(b) %*% A)
#   r <- sum(b*b)
#   y <- value(expr@args[[2]])
#   list(list(P/y), q/y, r/y)
# }
#
# .coeffs_power.QuadCoeffExtractor <- function(object, expr) {
#   if(expr@p == 1)
#     return(get_coeffs(object, expr@args[[1]]))
#   else if(expr@p == 2) {
#     coeffs <- .coeffs_affine(object, expr@args[[1]])
#     A <- coeffs[[2]]
#     b <- coeffs[[3]]
#
#     Ps <- lapply(1:nrow(A), function(i) { Matrix(t(A[i,]) %*% A[i,], sparse = TRUE) })
#     Q <- 2*Matrix(diag(b) %*% A, sparse = TRUE)
#     R <- b^2
#     return(list(Ps, Q, R))
#   } else
#     stop("Error while processing Power")
# }
#
# .coeffs_matrix_frac <- function(object, expr) {
#   coeffs <- .coeffs_affine(expr@args[[1]])
#   A <- coeffs[[2]]
#   b <- coeffs[[3]]
#   arg_size <- size(expr@args[[1]])
#   m <- arg_size[1]
#   n <- arg_size[2]
#   Pinv <- base::solve(value(expr@args[[2]]))
#
#   M <- sparseMatrix(i = c(), j = c(), dims = c(object@N, object@N))
#   Q <- sparseMatrix(i = c(), j = c(), dims = c(1, object@N))
#   R <- 0
#
#   for(i in seq(1, m*n, m)) {
#     A2 <- A[i:(i+m),]
#     b2 <- b[i:(i+m)]
#
#     M <- M + t(A2) %*% Pinv %*% A2
#     Q <- Q + 2*t(A2) %*% (Pinv %*% b2)
#     R <- R + sum(b2 * (Pinv %*% b2))
#   }
#   list(list(Matrix(M, sparse = TRUE)), Matrix(Q, sparse = TRUE), R)
# }
#
# .coeffs_affine_atom <- function(object, expr) {
#   sz <- prod(size(expr))
#   Ps <- lapply(1:sz, function(i) { sparseMatrix(i = c(), j = c(), dims = c(object@N, object@N)) })
#   Q <- sparseMatrix(i = c(), j = c(), dims = c(sz, object@N))
#   Parg <- NA
#   Qarg <- NA
#   Rarg <- NA
#
#   fake_args <- list()
#   offsets <- list()
#   offset <- 0
#   for(idx in 1:length(expr@args)) {
#     arg <- expr@args[[idx]]
#     if(is_constant(arg))
#       fake_args <- c(fake_args, list(create_const(value(arg), size(arg))))
#     else {
#       coeffs <- get_coeffs(object, arg)
#       if(is.na(Parg)) {
#         Parg <- coeffs[[1]]
#         Qarg <- coeffs[[2]]
#         Rarg <- coeffs[[3]]
#       } else {
#         Parg <- Parg + coeffs[[1]]
#         Qarg <- rbind(Qarg, q)
#         Rarg <- c(Rarg, r)
#       }
#       fake_args <- c(fake_args, list(create_var(size(arg), idx)))
#       offsets[idx] <- offset
#       offset <- offset + prod(size(arg))
#     }
#   }
#   graph <- graph_implementation(expr, fake_args, size(expr), get_data(expr))
#
#   # Get the matrix representation of the function
#   prob_mat <- get_problem_matrix(list(create_eq(graph[[1]])), offsets)
#   V <- prob_mat[[1]]
#   I <- prob_mat[[2]]
#   J <- prob_mat[[3]]
#   R <- as.vector(prob_mat[[4]])   # TODO: Check matrix is flattened correctly
#
#   # Return AX + b
#   for(idx in 1:length(V)) {
#     v <- V[idx]
#     i <- I[idx]
#     j <- J[idx]
#
#     Ps[[i]] <- Ps[[i]] + v*Parg[j]
#     Q[i,] <- Q[i,] + v*Qarg[j,]
#     R[i] <- R[i] + v*Rarg[j]
#   }
#   Ps <- lapply(Ps, function(P) { Matrix(P, sparse = TRUE) })
#   list(Ps, Matrix(Q, sparse = TRUE), R)
# }
