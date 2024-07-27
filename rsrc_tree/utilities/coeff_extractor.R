## CVXPY SOURCE: cvxpy/utilities/coef_extractor.py

# TODO: Find best format for sparse matrices.
.CoeffExtractor <- setClass("CoeffExtractor", representation(inverse_data = "InverseData", canon_backend = "character", id_map = "list", x_length = "numeric", var_dims = "list", param_dims = "list", param_to_size = "list", param_id_map = "list"),
                                              prototype(canon_backend = NA_character_, id_map = list(), x_length = NA_real_, var_dims = list(), param_dims = list(), param_to_size = list(), param_id_map = list()))

CoeffExtractor <- function(inverse_data, canon_backend = NA_character_) { .CoeffExtractor(inverse_data = inverse_data, canon_backend = canon_backend) }

setMethod("initialize", "CoeffExtractor", function(.Object, inverse_data, canon_backend = NA_character_, id_map = list(), x_length = NA_real_, var_dims = list(), param_dims = list(), param_to_size = list(), param_id_map = list()) {
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

      # TODO: All the sparse COO matrix calculations involving P, c_part, etc must be checked!
      c_part <- c[(var_offset + 1):(var_offset + var_size),]
      P <- value(quad_forms[[var_id_char]][[3]]@P)
      if(!any(is.na(P))) {
        # Convert to sparse matrix.
        P <- as(P, "TsparseMatrix")
        if(var_size == 1)
          c_part <- matrix(1, nrow = nrow(P), ncol = 1) * c_part
      } else
        P <- sparseMatrix(i = 1:var_size, j = 1:var_size, x = rep(1, var_size), repr = "T")

      # We multiply the columns of P by c_part by operating directly on the data.
      data <- as.matrix(P@x) * c_part[P@j]
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

          acc_data <- c(acc_data, data)
          acc_row <- c(acc_row, P@i)
          acc_col <- c(acc_col, P@j)

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

  # Extract quadratic matrices and vectors.
  num_params <- ncol(constant)
  P_list <- list()
  q_list <- list()
  P_height <- 0
  P_entries <- 0
  for(var_id_char in names(offsets)) {
    offset <- offsets[[var_id_char]]
    dim <- object@var_dims[[var_id_char]]
    size <- as.integer(prod(dim))

    if(var_id_char %in% names(coeffs) && "P" %in% names(coeffs[[var_id_char]])) {
      P <- coeffs[[var_id_char]]$P
      P_entries <- P_entries + length(P[[1]])
    } else
      P <- list(c(), list(c(), c()), c(size, size))

    if(var_id_char %in% names(coeffs) && "q" %in% coeffs[[var_id_char]])
      q <- coeffs[[var_id_char]]$q
    else
      q <- matrix(0, nrow = size, ncol = num_params)

    P_list <- c(P_list, list(P))
    q_list <- c(q_list, list(q))
    P_height <- P_height + size
  }

  if(P_height != object@x_length)
    stop("Resulting quadratic form does not have appropriate dimensions")

  # Conceptually we build a block diagonal matrix
  # out of all the Ps, then flatten the first two dimensions.
  # e.g. P1
  #        P2
  # We do this by extending each P with zero blocks above and below.
  gap_above <- bit64::integer64(0)
  acc_height <- bit64::integer64(0)
  rows <- as.integer64(rep(0, P_entries))
  cols <- as.integer64(rep(0, P_entries))
  vals <- rep(0, P_entries)
  entry_offset <- 0
  for(P in P_list) {
    ##Conceptually, the code is equivalent to
    ## above = np.zeros((gap_above, P.shape[1], num_params))
    ## below = np.zeros((gap_below, P.shape[1], num_params))
    ## padded_P = np.concatenate([above, P, below], axis=0)
    ## padded_P = np.reshape(padded_P, (P_height*P.shape[1], num_params),
    ##                       order='F')
    ## padded_P_list.append(padded_P)

    ## but done by constructing a COO matrix.

    padded_P_list.append(padded_P)
    P_vals <- P[[1]]
    P_rows <- P[[2]][[1]]
    P_cols <- P[[2]][[2]]
    P_dim <- P[[3]]

    if(!is.null(P_vals) && length(P_vals) > 0) {
      P_cols_ext <- big64::as.integer64(P_cols) * integer64(P_height)   # TODO: Check if this is elementwise multiplication of a rep call in Python.
      base_rows <- gap_above + acc_height + P_rows + P_cols_ext
      full_rows <- rep(base_rows, num_params)
      rows[(entry_offset + 1):(entry_offset + length(P_vals))] <- full_rows
      full_cols <- as.vector(matrix(seq(num_params), nrow = length(P_cols), ncol = num_params, byrow = TRUE))
      cols[(entry_offset + 1):(entry_offset + length(P_vals))] <- full_cols
      entry_offset <- entry_offset + length(P_vals)
    }

    gap_above <- gap_above + P_dim[1]
    acc_height <- acc_height + P_height*bit64::integer64(P_dim[2])
  }

  # Stitch together Ps and qs and constant.
  P <- sparseMatrix(i = rows, j = cols, x = vals, dims = c(acc_height, num_params))

  # Stack q with constant offset as last row.
  q <- do.call("rbind", q_list)
  q <- do.call("rbind", list(q, constant@A))
  q <- Matrix(q, sparse = TRUE)   # Store as CSR (row-wise) sparse matrix.
  return(list(P, q))
})
