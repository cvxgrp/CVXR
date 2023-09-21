# TODO: Find best format for sparse matrices.
.CoeffExtractor <- setClass("CoeffExtractor", representation(inverse_data = "InverseData", canon_backend = "character", id_map = "list", x_length = "integer", var_dims = "list", param_dims = "list", param_to_size = "list", param_id_map = "list"),
                                              prototype(canon_backend = NA_character_, id_map = list(), x_length = NA_integer_, var_dims = list(), param_dims = list(), param_to_size = list(), param_id_map = list()))

CoeffExtractor <- function(inverse_data, canon_backend = NA_character_) { .CoeffExtractor(inverse_data = inverse_data, canon_backend = canon_backend) }

setMethod("initialize", "CoeffExtractor", function(.Object, inverse_data, canon_backend = NA_character_, id_map = list(), x_length = NA_integer_, var_dims = list(), param_dims = list(), param_to_size = list(), param_id_map = list()) {
  .Object@inverse_data <- .Object@inverse_data
  .Object@canon_backend <- .Object@canon_backend
  .Object@id_map <- inverse_data@var_offsets
  .Object@x_length <- inverse_data@x_length
  .Object@var_dims <- inverse_data@var_dims
  .Object@param_dims <- inverse_data@param_dims
  .Object@param_to_size <- inverse_data@param_to_size
  .Object@param_id_map <- inverse_data@param_id_map
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
  expr_size <- size(expr)
  list(Matrix(nrow = expr_size, ncol = object@x_length), as.vector(value(expr)))
})

#' @describeIn CoeffExtractor Extract problem data tensor from an expression that is reducible to \eqn{A*x + b}. Applying the tensor to a flattened parameter vector and reshaping will recover \eqn{A} and \eqn{b} (see the helpers in canonInterface).
setMethod("affine", signature(object = "CoeffExtractor", expr = "list"), function(object, expr) {
  expr_list <- expr
  if(!all(sapply(expr_list, is_dpp)))
    stop("All expressions must be DPP")
  
  num_rows <- sum(sapply(expr_list, size))
  op_list <- lapply(expr_list, function(e) { canonical_form(e)[[1]] })
  canonInterface.get_problem_matrix(op_list, object@x_length, object@id_map, object@param_to_size, object@param_id_map, num_rows, object@canon_backend)
})

setMethod("affine", signature(object = "CoeffExtractor", expr = "Expression"), function(object, expr) {
  affine(object, list(expr))
})

setMethod("extract_quadratic_coeffs", "CoeffExtractor", function(object, affine_expr, quad_forms) {
  # Assumes quadratic forms all have variable arguments. Affine expressions can be anything.
  if(!is_dpp(affine_expr))
    stop("affine_expr must be DPP")
  
  # Here we take the problem objective, replace all the SymbolicQuadForm
  # atoms with variables of the same dimensions.
  # We then apply the canonInterface to reduce the "affine head"
  # of the expression tree to a coefficient vector c and constant offset d.
  # Because the expression is parameterized, we extend that to a matrix
  # [c1 c2 ...]
  # [d1 d2 ...]
  # where ci,di are the vector and constant for the ith parameter.
  tmp <- InverseData.get_var_offsets(variables(affine_expr))
  affine_id_map <- tmp[[1]]
  affine_offsets <- tmp[[2]]
  x_length <- tmp[[3]]
  affine_var_dims <- tmp[[4]]
  op_list <- list(canonical_form(affine_expr)[[1]])
  param_coeffs <- canonInterface.get_problem_matrix(op_list, x_length, affine_offsets, object@param_to_size, object@param_id_map, size(affine_expr), object@canon_backend)
  
  # Iterates over every entry of the parameters vector,
  # and obtains the Pi and qi for that entry i.
  # These are then combined into matrices [P1.flatten(), P2.flatten(), ...]
  # and [q1, q2, ...]
  constant <- param_coeffs[nrow(param_coeffs),]
  c <- param_coeffs[seq(nrow(param_coeffs) - 1),]@A
  
  # coeffs stores the P and q for each quad_form, as well as for true variable 
  # nodes in the objective.
  coeffs <- list()
  
  # The goal of this loop is to appropriately multiply the matrix P of each
  # quadratic term by the coefficients in param_coeffs. Later, we combine all 
  # the quadratic terms to form a single matrix P.
  for(var in variables(affine_expr)) {
    # quad_forms maps the ids of the SymbolicQuadForm atoms in the objective to
    # (modified parent node of quad form, argument index of quad form, quad 
    # form atom).
    var_id <- id(var)
    var_id_char <- as.character(vid)
    if(var_id_char %in% names(quad_forms)) {
      orig_id <- id(quad_forms[[var_id_char]][[3]]@args[[1]])
      var_offset <- affine_id_map[[var_id_char]][[1]]
      var_size <- affine_id_map[[var_id_char]][[2]]
      
      c_part <- c[(var_offset + 1):(var_offset + var_size),]
      P <- value(quad_forms[[var_id_char]][[3]]@P)
      if(!any(is.na(P))) {
        # Convert to sparse matrix.
        P <- as(P, "TsparseMatrix")
        if(var_size == 1)
          c_part <- matrix(1, nrow = nrow(P), ncol = 1) * c_part   # TODO: Check this multiplication is correct.
      } else
        P <- sparseMatrix(i = 1:var_size, j = 1:var_size, x = rep(1, var_size), repr = "T")
    
      # We multiply the columns of P by c_part by operating directly on the data.
      # TODO: Finish converting Python code.
      # data = P.data[:,None] * c_part[P.col]
      P_tup <- list(data, list(P@i, P@j), dim(P))
      
      # Conceptually similar to P = P[:,:,None] * c_part[None,:,:].
      orig_id_char <- as.character(orig_id)
      if(orig_id_char %in% names(coeffs)) {
        if("P" %in% names(coeffs[[orig_id_char]])) {
          # Concatenation becomes addition when constructing COO matrix because repeated indices are summed.
          # Conceptually equivalent to coeffs[[orig_id_char]]$P <- coeffs[[orig_id_char]]$P + P_tup.
          tmp <- coeffs[[orig_id_char]]$P
          acc_data <- tmp[[1]]
          acc_row <- tmp[[2]][[1]]
          acc_col <- tmp[[2]][[2]]
          
          # TODO: Finish converting Python code.
          # acc_data = np.concatenate([acc_data, data], axis = 0)
          # acc_row = np.concatenate([acc_row, P@i], axis = 0)
          # acc_col = np.concatenate([acc_col, P@j], axis = 0)
          
          P_tup <- list(acc_data, list(acc_row, acc_col), dim(P))
          coeffs[[orig_id_char]]$P <- P_tup
        } else
          coeffs[[orig_id_char]]$P <- P_tup
      } else {
        coeffs[[orig_id_char]] <- list()
        coeffs[[orig_id_char]]$P <- P_tup
        coeffs[[orig_id_char]]$q <- matrix(0, nrow = nrow(P), ncol = ncol(c))
      }
    } else {
      var_offset <- affine_id_map[[var_id_char]][[1]]
      var_size <- as.integer(prod(affine_var_dims[[var_id_char]]))
      if(var_id_char %in% names(coeffs))
        coeffs[[var_id_char]]$q <- coeffs[[var_id_char]]$q + c[(var_offset + 1):(var_offset + var_size),]
      else {
        coeffs[[var_id_char]] <- list()
        coeffs[[var_id_char]]$q <- c[(var_offset + 1):(var_offset + var_size),]
      }
    }
  }
  
  return(list(coeffs, constant))
  

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
        P <- sweep(P, MARGIN = 2, FUN = "*", c_part)
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