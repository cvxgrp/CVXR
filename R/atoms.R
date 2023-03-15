#'
#' The Atom class.
#'
#' This virtual class represents atomic expressions in CVXR.
#'
#' @name Atom-class
#' @aliases Atom
#' @rdname Atom-class
Atom <- setClass("Atom", representation(atom_args = "list", .dim = "NumORNULL"), prototype(atom_args = list(), .dim = NULL), contains = c("VIRTUAL", "Expression"))

setMethod("initialize", "Atom", function(.Object, ..., atom_args = list(), .dim = NULL, validate = TRUE) {
  # .Object@id <- ifelse(is.na(id), get_id(), id)
  .Object <- callNextMethod(.Object, ..., validate = FALSE)

  if(length(atom_args) == 0)
    stop("No arguments given to ", class(.Object), ".")
  .Object@args <- lapply(atom_args, as.Constant)
  validate_args(.Object)

  .Object@.dim <- dim_from_args(.Object)
  if(length(.Object@.dim) > 2)
    stop("Atoms must be at most 2D.")
  if(validate)
    validObject(.Object)
  .Object
})

setMethod("show", "Atom", function(object) {
  if(is.null(get_data(object)))
    data <- list()
  else
    data <- sapply(get_data(object), as.character)
  arg_names <- sapply(object@args, name)
  cat(class(object), "(", paste(c(arg_names, data), collapse = ", "), ")", sep = "")
})

#' @param x,object An \linkS4class{Atom} object.
#' @describeIn Atom Returns the string representtation of the function call
setMethod("name", "Atom", function(x) {
  if(is.null(get_data(x)))
    data <- list()
  else
    data <- sapply(get_data(x), as.character)
  arg_names <- sapply(x@args, name)
  paste(class(x), "(", paste(c(arg_names, data), collapse = ", "), ")", sep = "")
})

#' @describeIn Atom Raises an error if the arguments are invalid.
setMethod("validate_args", "Atom", function(object) {
  if(!allow_complex(object) && any(sapply(object@args, is_complex)))
    stop("Arguments to ", class(object), " cannot be complex.")
})

#' @rdname dim_from_args
setMethod("dim_from_args", "Atom", function(object) { stop("Unimplemented") })

#' @describeIn Atom The \code{c(row, col)} dimensions of the atom.
setMethod("dim", "Atom", function(x) { x@.dim })

#' @describeIn Atom The number of rows in the atom.
setMethod("nrow", "Atom", function(x) { dim(x)[1] })

#' @describeIn Atom The number of columns in the atom.
setMethod("ncol", "Atom", function(x) { dim(x)[2] })

#' @describeIn Atom Does the atom handle complex numbers?
setMethod("allow_complex", "Atom", function(object) { FALSE })

#' @rdname sign_from_args
setMethod("sign_from_args", "Atom", function(object) { stop("Unimplemented") })

#' @describeIn Atom A logical value indicating whether the atom is nonnegative.
setMethod("is_nonneg", "Atom", function(object) { sign_from_args(object)[1] })

#' @describeIn Atom A logical value indicating whether the atom is nonpositive.
setMethod("is_nonpos", "Atom", function(object) { sign_from_args(object)[2] })

#' @describeIn Atom A logical value indicating whether the atom is imaginary.
setMethod("is_imag", "Atom", function(object) { FALSE })

#' @describeIn Atom A logical value indicating whether the atom is complex valued.
setMethod("is_complex", "Atom", function(object) { FALSE })

#' @rdname curvature-atom
setMethod("is_atom_convex", "Atom", function(object) { stop("Unimplemented") })

#' @rdname curvature-atom
setMethod("is_atom_concave", "Atom", function(object) { stop("Unimplemented") })

#' @rdname curvature-atom
setMethod("is_atom_affine", "Atom", function(object) { is_atom_concave(object) && is_atom_convex(object) })

#' @rdname curvature-atom
setMethod("is_atom_log_log_convex", "Atom", function(object) { FALSE })

#' @rdname curvature-atom
setMethod("is_atom_log_log_concave", "Atom", function(object) { FALSE })

#' @rdname curvature-atom
setMethod("is_atom_log_log_affine", "Atom", function(object) { is_atom_log_log_concave(object) && is_atom_log_log_convex(object) })

#' @rdname curvature-atom
setMethod("is_atom_quasiconvex", "Atom", function(object) { is_atom_convex(object) })

#' @rdname curvature-atom
setMethod("is_atom_quasiconcave", "Atom", function(object) { is_atom_concave(object) })

#' @rdname curvature-comp
setMethod("is_incr", "Atom", function(object, idx) { stop("Unimplemented") })

#' @rdname curvature-comp
setMethod("is_decr", "Atom", function(object, idx) { stop("Unimplemented") })

#' @describeIn Atom A logical value indicating whether the atom is convex.
setMethod("is_convex", "Atom", function(object) {
  # Applies DCP composition rule
  if(is_constant(object))
    return(TRUE)
  else if(is_atom_convex(object)) {
    idx <- 1
    for(arg in object@args) {
      if(!(is_affine(arg) || (is_convex(arg) && is_incr(object, idx)) || (is_concave(arg) && is_decr(object, idx))))
        return(FALSE)
      idx <- idx + 1
    }
    return(TRUE)
  } else
    return(FALSE)
})

#' @describeIn Atom A logical value indicating whether the atom is concave.
setMethod("is_concave", "Atom", function(object) {
  # Applies DCP composition rule
  if(is_constant(object))
    return(TRUE)
  else if(is_atom_concave(object)) {
    idx <- 1
    for(arg in object@args) {
      if(!(is_affine(arg) || (is_concave(arg) && is_incr(object, idx)) || (is_convex(arg) && is_decr(object, idx))))
        return(FALSE)
      idx <- idx + 1
    }
    return(TRUE)
  } else
    return(FALSE)
})

#' @describeIn Atom A logical value indicating whether the atom is a disciplined parameterized expression.
setMethod("is_dpp", "Atom", function(object, context = "dcp") {
  if(tolower(context) == "dcp")
    return(is_dcp(object, dpp = TRUE))
  else if(tolower(context) == "dgp")
    return(is_dgp(object, dpp = TRUE))
  else:
    stop("Unsupported context ", context)
})

#' @describeIn Atom A logical value indicating whether the atom is log-log convex.
setMethod("is_log_log_convex", "Atom", function(object) {
  # Verifies DGP composition rule.
  if(is_log_log_constant(object))
    return(TRUE)
  else if(is_atom_log_log_convex(object)) {
    idx <- 1
    for(arg in object@args) {
      if(!(is_log_log_affine(arg) || (is_log_log_convex(arg) && is_incr(object, idx)) || (is_log_log_concave(arg) && is_decr(object, idx))))
        return(FALSE)
      idx <- idx + 1
    }
    return(TRUE)
  } else
    return(FALSE)
})

#' @describeIn Atom A logical value indicating whether the atom is log-log concave.
setMethod("is_log_log_concave", "Atom", function(object) {
  # Verifies DGP composition rule.
  if(is_log_log_constant(object))
    return(TRUE)
  else if(is_atom_log_log_concave(object)) {
    idx <- 1
    for(arg in object@args) {
      if(!(is_log_log_affine(arg) || (is_log_log_concave(arg) && is_incr(object, idx)) || (is_log_log_convex(arg) && is_decr(object, idx))))
        return(FALSE)
      idx <- idx + 1
    }
    return(TRUE)
  } else
    return(FALSE)
})

setMethod(".non_const_idx", "Atom", function(object) {
  return(which(sapply(object@args, is_constant)))
})

setMethod(".is_real", "Atom", function(object) {
  # Returns TRUE if this atom is a real function:
  #    The atom must have exactly one argument that is not a constant.
  #    The argument must be a scalar.
  #    The output must be a scalar.
  non_const <- .non_const_idx(object)
  return(is_scalar(object) && len(non_const) == 1 && is_scalar(object@args[[non_const[1]]]))
})

#' @describeIn Atom A logical value indicating whether the atom is quasiconvex.
setMethod("is_quasiconvex", "Atom", function(object) {
  if(is_convex(object))
    return(TRUE)
  if(class(object) == "MaxElemwise" || class(object) == "MaxEntries")
    return(all(sapply(object@args, is_quasiconvex)))
  non_const <- .non_const_idx(object)
  if(.is_real(object) && is_incr(object, non_const[[1]]))
    return(is_quasiconvex(object@args[[non_const[[1]]]]))
  if(.is_real(object) && is_decr(object, non_const[[1]]))
    return(is_quasiconcave(object@args[[non_const[[1]]]]))
  if(is_atom_quasiconvex(object)) {
    idx <- 1
	for(arg in object@args) {
	  if(!(is_affine(arg) || (is_convex(arg) && is_incr(object, idx)) || (is_concave(arg) && is_decr(object, idx))))
	    return(FALSE)
	  idx <- idx + 1
    }
	return(TRUE)
  }
  return(FALSE)
})

#' @describeIn Atom A logical value indicating whether the atom is quasiconcave.
setMethod("is_quasiconcave", "Atom", function(object) {
  if(is_concave(object))
    return(TRUE)
  if(class(object) == "MinElemwise" || class(object) == "MinEntries")
    return(all(sapply(object@args, is_quasiconcave)))
  non_const <- .non_const_idx(object)
  if(.is_real(object) && is_incr(object, non_const[[1]]))
    return(is_quasiconcave(object@args[[non_const[[1]]]]))
  if(.is_real(object) && is_decr(object, non_const[[1]]))
    return(is_quasiconvex(object@args[[non_const[[1]]]]))
  if(is_atom_quasiconcave(object)) {
    idx <- 1
	for(arg in object@args) {
	  if(!(is_affine(arg) || (is_concave(arg) && is_incr(object, idx)) || (is_convex(arg) && is_decr(object, idx))))
	    return(FALSE)
	  idx <- idx + 1
    }
	return(TRUE)
  }
  return(FALSE)
})

#' @describeIn Atom Represent the atom as an affine objective and conic constraints.
setMethod("canonicalize", "Atom", function(object) {
  # Constant atoms are treated as a leaf.
  if(is_constant(object) && !is.na(parameters(object)) && length(parameters(object)) > 0)
      # Non-parameterized expressions are evaluated immediately.
      return(canonical_form(Constant(value(object))))
  else {
    arg_objs <- list()
    constraints <- list()
    for(arg in object@args) {
      canon <- canonical_form(arg)
      arg_objs[[length(arg_objs) + 1]] <- canon[[1]]
      constraints <- c(constraints, canon[[2]])
    }
    # Special info required by the graph implementation.
    data <- get_data(object)
    graph <- graph_implementation(object, arg_objs, dim(object), data)
    return(list(graph[[1]], c(constraints, graph[[2]])))
  }
})

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector with two elements representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Atom The graph implementation of the atom.
setMethod("graph_implementation", "Atom", function(object, arg_objs, dim, data = NA_real_) { stop("Unimplemented") })

# .value_impl.Atom <- function(object) {
#' @describeIn Atom Returns the value of each of the componets in an Atom. Returns an empty matrix if it's an empty atom
setMethod("value_impl", "Atom", function(object) {
  obj_dim <- dim(object)
  # dims with 0's dropped in presolve.
  if(0 %in% obj_dim)
    result <- matrix(nrow = 0, ncol = 0)
  else {
    arg_values <- list()
    for(arg in object@args) {
      # An argument without a value makes all higher level values NA.
      # But if the atom is constant with non-constant arguments, it doesn't depend on its arguments, so it isn't NA.
      arg_val <- value_impl(arg)
      if(any(is.na(arg_val)) && !is_constant(object))
        return(NA_real_)
      else
        arg_values <- c(arg_values, list(arg_val))
    }
    result <- to_numeric(object, arg_values)
  }
  return(result)
})

#' @describeIn Atom Returns the value of the atom.
setMethod("value", "Atom", function(object) {
  if(any(sapply(parameters(object), function(p) { is.na(value(p)) })))
    return(NA_real_)
  return(value_impl(object))
})

#' @describeIn Atom The (sub/super)-gradient of the atom with respect to each variable.
setMethod("grad", "Atom", function(object) {
  # Short-circuit to all zeros if known to be constant
  if(is_constant(object))
    return(constant_grad(object))

  # Returns NA if variable values are not supplied
  arg_values <- list()
  for(arg in object@args) {
    arg_val <- value(arg)
    if(any(is.na(arg_val)))
      return(error_grad(object))
    else
      arg_values <- c(arg_values, list(arg_val))
  }

  # A list of gradients wrt arguments
  grad_self <- .grad(object, arg_values)

  # The chain rule
  result <- list()
  idx <- 1
  for(arg in object@args) {
    # A dictionary of gradients wrt variables.
    # Partial argument / partial x.
    grad_arg <- grad(arg)
    for(key in names(grad_arg)) {
      # None indicates gradient is not defined.
      if(any(is.na( as.vector(grad_arg[[key]]) )) || any(is.na( as.vector(grad_self[[idx]]) )))
        result[[key]] <- NA_real_
      else {
        D <- grad_arg[[key]] %*% grad_self[[idx]]
        # Convert 1x1 matrices to scalars.
        if((is.matrix(D) || is(D, "Matrix")) && all(dim(D) == c(1,1)))
          D <- D[1,1]

        if(key %in% names(result))
          result[[key]] <- result[[key]] + D
        else
          result[[key]] <- D
      }
    }
    idx <- idx + 1
  }
  return(result)
})

setMethod(".grad", "Atom", function(object, values) { stop("Unimplemented") })

#' @describeIn Atom A list of constraints describing the closure of the region where the expression is finite.
setMethod("domain", "Atom", function(object) {
  cons <- list()
  for(arg in object@args) {
    for(con in domain(arg))
      cons <- c(cons, con)
  }
  c(.domain(object), cons)
})

setMethod(".domain", "Atom", function(object) { list() })

#' @describeIn Atom Returns a list of the atom types present amongst this atom's arguments
setMethod("atoms", "Atom", function(object) {
  atom_list <- list()
  for(arg in object@args)
    atom_list <- c(atom_list, atoms(arg))
  atom_list <- c(atom_list, list(class(object)))
  return(unique(atom_list))
})

#' @param value A numeric value
#' @describeIn AxisAtom Gives the (sub/super)gradient of the atom w.r.t. each column variable
setMethod(".column_grad", "AxisAtom", function(object, value) { stop("Unimplemented") })

#'
#' The AxisAtom class.
#'
#' This virtual class represents atomic expressions that can be applied along an axis in CVXR.
#'
#' @slot expr A numeric element, data.frame, matrix, vector, or Expression.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @slot keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @name AxisAtom-class
#' @aliases AxisAtom
#' @rdname AxisAtom-class
AxisAtom <- setClass("AxisAtom", representation(expr = "ConstValORExpr", axis = "ANY", keepdims = "logical"),
                                 prototype(axis = NA_real_, keepdims = FALSE), contains = c("VIRTUAL", "Atom"))

setMethod("initialize", "AxisAtom", function(.Object, ..., expr, axis = NA_real_, keepdims = FALSE) {
  .Object@expr <- expr
  .Object@axis <- axis
  .Object@keepdims <- keepdims
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object An \linkS4class{Atom} object.
#' @describeIn AxisAtom The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "AxisAtom", function(object) {
  # TODO: Revisit this when we properly handle dimensions of scalars (NULL) and 1-D vectors (length only).
  arg_dim <- dim(object@args[[1]])
  if(object@keepdims && is.na(object@axis))   # Copy scalar to maintain original dimensions.
    arg_dim <- rep(1, length(arg_dim))
  else if(object@keepdims && !is.na(object@axis)) {   # Collapse dimensions NOT in axis to 1.
    collapse <- setdiff(1:length(arg_dim), object@axis)
    arg_dim[collapse] <- 1
  } else if(!object@keepdims && is.na(object@axis))   # Return a scalar.
    # arg_dim <- NULL   # TODO: Should this be NA instead?
    arg_dim <- rep(1, length(arg_dim))
  else {   # Drop dimensions NOT in axis and collapse atom.
    # arg_dim <- arg_dim[object@axis]
    collapse <- setdiff(1:length(arg_dim), object@axis)
    arg_dim <- c(arg_dim[object@axis], rep(1, length(collapse)))
  }
  return(arg_dim)
})

#' @describeIn AxisAtom A list containing \code{axis} and \code{keepdims}.
setMethod("get_data", "AxisAtom", function(object) { list(object@axis, object@keepdims) })

#' @describeIn AxisAtom Check that the new dimensions have the same number of entries as the old.
setMethod("validate_args", "AxisAtom", function(object) {
  if(!is.na(object@axis) && any(object@axis > ndim(object@args[[1]]) || object@axis <= 0))
    stop("Invalid argument for axis. Must be an integer between 1 and ", ndim(object@args[[1]]))
  callNextMethod()
})

#' @param values A list of numeric values for the arguments
#' @describeIn AxisAtom Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".axis_grad", "AxisAtom", function(object, values) {
  if(is.na(object@axis) || ndim(object@args[[1]]) < 2) {
    value <- matrix(values[[1]], nrow = size(object@args[[1]]), ncol = 1)
    D <- .column_grad(object, value)
    if(is(D, "Matrix") || !any(is.na(D)))
      D <- Matrix(D, sparse = TRUE)
  } else {
    m <- nrow(object@args[[1]])
    n <- ncol(object@args[[1]])
    if(object@axis == 2) {   # Function apply to each column
      D <- sparseMatrix(i = c(), j = c(), dims = c(m*n, n))
      for(i in 1:n) {
        value <- values[[1]][,i]
        d <- t(.column_grad(object, value))
        if(any(is.na(as.vector(d))))
          return(list(NA_real_))
        row <- seq((i-1)*n+1, (i-1)*n+m, length.out = m)
        col <- rep(1,m) * i
        D <- D + sparseMatrix(i = row, j = col, x = as.vector(d), dims = c(m*n, n))
      }
    } else {   # Function apply to each row
      values <- t(values[[1]])
      D <- sparseMatrix(i = c(), j = c(), dims = c(m*n, m))
      for(i in 1:m) {
        value <- values[,i]
        d <- t(.column_grad(object, value))
        if(any(is.na(as.vector(d))))
          return(list(NA_real_))
        row <- seq(i, i+(n-1)*m, length.out = n)
        col <- rep(1,n)*i
        D <- D + sparseMatrix(i = row, j = col, x = as.vector(d), dims = c(m*n, m))
      }
    }
  }
  list(D)
})

#'
#' The ConditionNumber class.
#'
#' This class represents the condition number of a positive definite matrix \eqn{A}, i.e.,
#'
#' \deqn{\lambda_{\max}(A) / \lambda_{\min}(A)}
#'
#' @slot expr An \linkS4class{Expression} representing a positive definite matrix.
#' @name ConditionNumber-class
#' @aliases ConditionNumber
#' @rdname ConditionNumber-class
.ConditionNumber <- setClass("ConditionNumber", representation(A = "ConstValORExpr"), contains = "Atom")

#' @param A An \linkS4class{Expression} or numeric matrix.
#' @rdname ConditionNumber-class
ConditionNumber <- function(A) { .ConditionNumber(A = A) }

setMethod("initialize", "ConditionNumber", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., atom_args = list(.Object@A))
})

#' @param object A \linkS4class{ConditionNumber} object.
#' @param values A list of arguments to the atom.
#' @describeIn ConditionNumber The condition number of the matrix.
setMethod("to_numeric", "ConditionNumber", function(object, values) {
  # eigen_vec <- base::eigen(values[[1]], only.values = TRUE)$values
  # max_eigen <- max(eigen_vec)
  # min_eigen <- min(eigen_vec)
  # return(max_eigen / min_eigen)
  base::kappa(values[[1]], exact = TRUE)
})

#' @describeIn ConditionNumber Returns constraints describing the domain of the node.
setMethod(".domain", "ConditionNumber", function(object) { list(Conj(t(object@args[[1]])) == object@args[[1]], object@args[[1]] %>>% 0) })

#' @param values A list of numeric values for the arguments
#' @describeIn ConditionNumber Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "ConditionNumber", function(object, values) { stop("Unimplemented") })

#' @describeIn ConditionNumber Check that the matrix is square.
setMethod("validate_args", "ConditionNumber", function(object) {
  arg_dim <- dim(object@args[[1]])
  if(length(arg_dim) != 2 || arg_dim[1] != arg_dim[2])
    stop("The argument to ConditionNumber must be a square matrix")
})

#' @describeIn ConditionNumber The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "ConditionNumber", function(object) { c(1,1) })

#' @describeIn ConditionNumber The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "ConditionNumber", function(object) { c(TRUE, FALSE) })

#' @describeIn ConditionNumber Is the atom convex?
setMethod("is_atom_convex", "ConditionNumber", function(object) { FALSE })

#' @describeIn ConditionNumber Is the atom concave?
setMethod("is_atom_concave", "ConditionNumber", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn ConditionNumber Is the atom weakly increasing in the index?
setMethod("is_incr", "ConditionNumber", function(object, idx) { FALSE })

#' @describeIn ConditionNumber Is the atom weakly decreasing in the index?
setMethod("is_decr", "ConditionNumber", function(object, idx) { FALSE })

#'
#' The CumMax class.
#'
#' This class represents the cumulative maximum of an expression.
#'
#' @slot expr An \linkS4class{Expression}.
#' @slot axis A numeric vector indicating the axes along which to apply the function. For a 2D matrix, \code{1} indicates rows, \code{2} indicates columns, and \code{c(1,2)} indicates rows and columns.
#' @name CumMax-class
#' @aliases CumMax
#' @rdname CumMax-class
.CumMax <- setClass("CumMax", prototype = prototype(axis = 2), contains = "AxisAtom")

#' @param expr An \linkS4class{Expression}.
#' @param axis A numeric vector indicating the axes along which to apply the function. For a 2D matrix, \code{1} indicates rows, \code{2} indicates columns, and \code{c(1,2)} indicates rows and columns.
#' @rdname CumMax-class
CumMax <- function(expr, axis = 2) { .CumMax(expr = expr, axis = axis) }

#' @param object A \linkS4class{CumMax} object.
#' @param values A list of arguments to the atom.
#' @describeIn CumMax The cumulative maximum along the axis.
setMethod("to_numeric", "CumMax", function(object, values) {
  # apply(values[[1]], object@axis, base::cummax)
  if(object@axis == 1)
    do.call(rbind, lapply(seq_len(nrow(values[[1]])), function(i) { base::cummax(values[[1]][i,]) }))
  else if(object@axis == 2)
    do.call(cbind, lapply(seq_len(ncol(values[[1]])), function(j) { base::cummax(values[[1]][,j]) }))
  else
    base::cummax(values[[1]])
})

#' @param values A list of numeric values for the arguments
#' @describeIn CumMax Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "CumMax", function(object, values) { .axis_grad(object, values) })

#' @param value A numeric value.
#' @describeIn CumMax Gives the (sub/super)gradient of the atom w.r.t. each column variable
setMethod(".column_grad", "CumMax", function(object, value) {
  # Grad: 1 for a largest index.
  value <- as.vector(value)
  maxes <- base::cummax(value)
  D <- matrix(0, nrow = length(value), ncol = 1)
  D[1,1] <- 1
  if(length(value) > 1)
    D[2:nrow(D),] <- maxes[2:length(maxes)] > maxes[1:(length(maxes)-1)]
  return(D)
})

#' @describeIn CumMax The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "CumMax", function(object) { dim(object@args[[1]]) })

#' @describeIn CumMax The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "CumMax", function(object) { c(is_nonneg(object@args[[1]]), is_nonpos(object@args[[1]])) })

#' @describeIn CumMax Returns the axis along which the cumulative max is taken.
setMethod("get_data", "CumMax", function(object) { list(object@axis) })

#' @describeIn CumMax Is the atom convex?
setMethod("is_atom_convex", "CumMax", function(object) { TRUE })

#' @describeIn CumMax Is the atom concave?
setMethod("is_atom_concave", "CumMax", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn CumMax Is the atom weakly increasing in the index?
setMethod("is_incr", "CumMax", function(object, idx) { TRUE })

#' @describeIn CumMax Is the atom weakly decreasing in the index?
setMethod("is_decr", "CumMax", function(object, idx) { FALSE })

#'
#' The DistRatio class.
#'
#' This class represents the ratio of the \eqn{l_2}-norm distance from \eqn{x} to two points \eqn{a} and \eqn{b}, i.e.,
#'
#' \deqn{||x - a||_2 / ||x - b||_2},
#' 
#' where \eqn{a} and \eqn{b} are constants with \eqn{norm(x - a)_2 \leq norm(x - b)}.
#'
#' @slot x An \linkS4class{Expression}.
#' @name DistRatio-class
#' @aliases DistRatio
#' @rdname DistRatio-class
.DistRatio <- setClass("DistRatio", representation(x = "ConstValORExpr", a = "numeric", b = "numeric"), 
                       validity = function(object) {
                         if(!is_constant(object@args[[2]]))
                           stop("[DistRatio: a] The argument a must be a constant.")
                         if(!is_constant(object@args[[3]]))
                           stop("[DistRatio: b] The argument b must be a constant.")
                         return(TRUE)
                       }, contains = "Atom")

#' @param x An \linkS4class{Expression}.
#' @param a A numeric value.
#' @param b A numeric value.
#' @rdname DistRatio-class
DistRatio <- function(x, a, b) { .DistRatio(x = x, a = a, b = b) }

setMethod("initialize", "DistRatio", function(.Object, ..., x, a, b) {
  .Object@x <- x
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@x, a, b))
  .Object@a <- value(.Object@args[[2]])
  .Object@b <- value(.Object@args[[3]])
  .Object
})

#' @param object A \linkS4class{DistRatio} object.
#' @param values A list of arguments to the atom.
#' @describeIn DistRatio The distance ratio of the expression.
setMethod("to_numeric", "DistRatio", function(object, values) {
  sqrt(sum((values[[1]] - object@a)^2)) / sqrt(sum((values[[1]] - object@b)^2))
})

#' @describeIn DistRatio The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "DistRatio", function(object) { c(1,1) })

#' @describeIn DistRatio The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "DistRatio", function(object) { c(TRUE, FALSE) })

#' @describeIn DistRatio Is the atom convex?
setMethod("is_atom_convex", "DistRatio", function(object) { FALSE })

#' @describeIn DistRatio Is the atom concave?
setMethod("is_atom_concave", "DistRatio", function(object) { FALSE })

#' @describeIn DistRatio Is the atom quasiconvex?
setMethod("is_atom_quasiconvex", "DistRatio", function(object) { TRUE })

#' @describeIn DistRatio Is the atom quasiconcave?
setMethod("is_atom_quasiconcave", "DistRatio", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn DistRatio Is the atom weakly increasing in the index?
setMethod("is_incr", "DistRatio", function(object, idx) { FALSE })

#' @describeIn DistRatio Is the atom weakly decreasing in the index?
setMethod("is_decr", "DistRatio", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn DistRatio Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "DistRatio", function(object, values) { return(NA_real_) })

#'
#' The DotSort class.
#'
#' This class represents the value
#' 
#' \deqn{\langle sort\left(vec(X)\right), sort\left(vec(W)\right) \rangle},
#' 
#' where \eqn{X} is an expression and \eqn{W} is a constant.
#' 
#' Both arguments are flattened, i.e., we define \eqn{x = vec(X)} and \eqn{w = vec(W)}.
#' If the length of \eqn{w} is less than the length of \eqn{x}, it is conceptually padded with zeroes.
#' When the length of \eqn{w} is larger than the length of \eqn{x}, an exception is raised.
#' 
#' DotSort is a generalization of SumLargest and SumSmallest:
#' SumLargest(X, 3) is equivalent to DotSort(X, c(1,1,1))
#' SumLargest(X, 3.5) is equivalent to DotSort(X, c(1,1,1,0.5))
#' SumSmallest(X, 3) is equivalent to -DotSort(X, c(-1,-1,-1))
#' 
#' When the constant argument is not a boolean vector, DotSort can be considered as a weighted sum 
#' of \eqn{x}, where the largest weight is assigned to the largest entry in \eqn{x}, etc.
#'
.DotSort <- setClass("DotSort", representation(X = "ConstValORExpr", W = "ConstVal"), contains = "Atom")

#' @param X An \linkS4class{Expression}.
#' @param W A numeric matrix.
#' @rdname DotSort-class
DotSort <- function(X, W) { .DotSort(X = X, W = W) }

setMethod("initialize", "DotSort", function(.Object, ..., X, W) {
  .Object@X <- X
  .Object@W <- W
  callNextMethod(.Object, ..., atom_args = list(.Object@X, .Object@W))
})

#' @describeIn DotSort Check that the arguments are valid.
setMethod("validate_args", "DotSort", function(object) {
  if(!is_constant(object@args[[2]]))
    stop("W must be constant")
  if(size(object@args[[1]]) < size(object@args[[2]]))
    stop("The size of W must be less than or equal to the size of X")
  callNextMethod()
})

#' @param object A \linkS4class{DotSort} object.
#' @param values A list of arguments to the atom.
#' @describeIn DotSort The inner product of the sorted values of vec(X) and the sorted (and potentially padded) values of vec(W).
setMethod("to_numeric", "DotSort", function(object, values) {
  args <- DotSort.get_args_from_values(values)
  x <- args[[1]]
  w_padded <- args[[2]]
  return(as.numeric(sort(x) %*% sort(w_padded)))
})

#' @param values A list of numeric values for the arguments
#' @describeIn DotSort Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "DotSort", function(object, values) {
  args <- DotSort.get_args_from_values(values)
  x <- args[[1]]
  w_padded <- args[[2]]
  n <- length(x)
  sorted_w <- sort(w_padded)
  return(list(sparseMatrix(i = indices, j = rep(1, n), x = sorted_w, dims = c(n, 1))))
})

#' @describeIn DotSort The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "DotSort", function(object) { c(1,1) })

#' @describeIn DotSort The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "DotSort", function(object) {
  # Same as argument.
  x_pos <- is_nonneg(object@args[[1]])
  x_neg <- is_nonpos(object@args[[1]])
  
  w_pos <- is_nonneg(object@args[[2]])
  w_neg <- is_nonpos(object@args[[2]])
  
  is_positive <- (x_pos && w_pos) || (x_neg && w_neg)
  is_negative <- (x_neg && w_pos) || (x_pos && w_neg)
  
  return(c(is_positive, is_negative))
})

#' @describeIn DotSort Is the atom convex?
setMethod("is_atom_convex", "DotSort", function(object) {
  if(dpp_scope_active()) {   # TODO: Figure out how to save DPP status as global parameter.
    # DotSort is convex under DPP if W is parameter affine.
    X <- object@args[[1]]
    W <- object@args[[2]]
    return(is_constant(X) || is_param_affine(W))
  } else
    return(TRUE)
})

#' @describeIn DotSort Is the atom concave?
setMethod("is_atom_concave", "DotSort", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn DotSort Is the atom weakly increasing in the index?
setMethod("is_incr", "DotSort", function(object, idx) { is_nonneg(object@args[[2]]) })

#' @describeIn DotSort Is the atom weakly decreasing in the index?
setMethod("is_decr", "DotSort", function(object, idx) { is_nonpos(object@args[[2]]) })

#' @describeIn DotSort Empty list because W is stored as an argument.
setMethod("get_data", "DotSort", function(object) { list() })

DotSort.get_args_from_values <- function(values) {
  x <- as.numeric(t(values[[1]]))
  w <- as.numeric(t(values[[2]]))
  
  w_padded <- rep(0, length(x))
  w_padded[1:length(w)] <- w
  return(list(x, w_padded))
}

#'
#' The EyeMinusInv class.
#'
#' This class represents the unity resolvent of an elementwise positive matrix \eqn{X}, i.e., \eqn{(I - X)^{-1}},
#' and it enforces the constraint that the spectral radius of \eqn{X} is at most 1.
#' This atom is log-log convex.
#'
#' @slot X An \linkS4class{Expression} or numeric matrix.
#' @name EyeMinusInv-class
#' @aliases EyeMinusInv
#' @rdname EyeMinusInv-class
.EyeMinusInv <- setClass("EyeMinusInv", representation(X = "ConstValORExpr"),
                         validity = function(object) {
                           if(length(dim(object@X)) != 2 || nrow(object@X) != ncol(object@X))
                             stop("[EyeMinusInv: X] The argument X must be a square matrix.")
                           return(TRUE)
                          }, contains = "Atom")

#' @param X An \linkS4class{Expression} or numeric matrix.
#' @rdname EyeMinusInv-class
EyeMinusInv <- function(X) { .EyeMinusInv(X = X) }

setMethod("initialize", "EyeMinusInv", function(.Object, ..., X) {
  .Object@X <- X
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@X))
  .Object@args[[1]] <- X
  .Object
})

#' @param object,x An \linkS4class{EyeMinusInv} object.
#' @param values A list of arguments to the atom.
#' @describeIn EyeMinusInv The unity resolvent of the matrix.
setMethod("to_numeric", "EyeMinusInv", function(object, values) {
  base::solve(diag(nrow(object@args[[1]])) - values[[1]])
})

#' @describeIn EyeMinusInv The name and arguments of the atom.
setMethod("name", "EyeMinusInv", function(x) { paste(class(x), x@args[[1]]) })

#' @describeIn EyeMinusInv The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "EyeMinusInv", function(object) { dim(object@args[[1]]) })

#' @describeIn EyeMinusInv The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "EyeMinusInv", function(object) { c(TRUE, FALSE) })

#' @describeIn EyeMinusInv Is the atom convex?
setMethod("is_atom_convex", "EyeMinusInv", function(object) { FALSE })

#' @describeIn EyeMinusInv Is the atom concave?
setMethod("is_atom_concave", "EyeMinusInv", function(object) { FALSE })

#' @describeIn EyeMinusInv Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "EyeMinusInv", function(object) { TRUE })

#' @describeIn EyeMinusInv Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "EyeMinusInv", function(object) { FALSE })

# TODO: Figure out monotonicity.
#' @param idx An index into the atom.
#' @describeIn EyeMinusInv Is the atom weakly increasing in the index?
setMethod("is_incr", "EyeMinusInv", function(object, idx) { FALSE })

#' @describeIn EyeMinusInv Is the atom weakly decreasing in the index?
setMethod("is_decr", "EyeMinusInv", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn EyeMinusInv Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "EyeMinusInv", function(object, values) { NA_real_ })

# The resolvent of a positive matrix, (sI - X)^(-1).
# For an elementwise positive matrix X and a positive scalar s, this atom computes
# (sI - X)^(-1), and it enforces the constraint that the spectral radius of X/s is
# at most 1.
# This atom is log-log convex.
Resolvent <- function(X, s) {
  1.0 / (s * EyeMinusInv(X / s))
}

#'
#' The GenLambdaMax class.
#'
#' This class represents the maximum generalized eigenvalue  \eqn{\lambda_{\max}(A, B)},
#' where \eqn{A} is a symmetric matrix and \eqn{B} is a positive semidefinite matrix.
#'
#' @slot A An \linkS4class{Expression} representing a symmetric matrix.
#' @slot B An \linkS4class{Expression} representing a positive semidefinite matrix.
#' @name GenLambdaMax-class
#' @aliases GenLambdaMax
#' @rdname GenLambdaMax-class
.GenLambdaMax <- setClass("GenLambdaMax", representation(A = "ConstValORExpr", B = "ConstValORExpr"), contains = "Atom")

#' @param A An \linkS4class{Expression} or numeric matrix.
#' @param B An \linkS4class{Expression} or numeric matrix.
#' @rdname GenLambdaMax-class
GenLambdaMax <- function(A, B) { .GenLambdaMax(A = A, B = B) }

setMethod("initialize", "GenLambdaMax", function(.Object, ..., A, B) {
  .Object@A <- A
  .Object@B <- B
  callNextMethod(.Object, ..., atom_args = list(.Object@A, .Object@B))
})

#' @param object A \linkS4class{GenLambdaMax} object.
#' @param values A list of arguments to the atom.
#' @describeIn GenLambdaMax The largest generalized eigenvalue corresponding to the matrices.
setMethod("to_numeric", "GenLambdaMax", function(object, values) {
  eigen_res <- geigen(values[[1]], values[[2]], symmetric = TRUE, only.values = TRUE)
  base::max(eigen_res$values)
})

#' @describeIn GenLambdaMax Returns constraints describing the domain of the node.
setMethod(".domain", "GenLambdaMax", function(object) { list(Conj(t(object@args[[1]])) == object@args[[1]], Conj(t(object@args[[2]])) == object@args[[2]], object@args[[2]] %>>% 0) })

#' @param values A list of numeric values for the arguments
#' @describeIn GenLambdaMax Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "GenLambdaMax", function(object, values) { stop("Unimplemented") })

#' @describeIn GenLambdaMax Check that the matrices are square and of the same dimension.
setMethod("validate_args", "GenLambdaMax", function(object) {
  A_dim <- dim(object@args[[1]])
  B_dim <- dim(object@args[[2]])
  if(length(A_dim) != 2 || A_dim[1] != A_dim[2] || B_dim[1] != B_dim[2] || !all(A_dim == B_dim))
    stop("The arguments to GenLambdaMax must be square and have the same dimensions")
})

#' @describeIn GenLambdaMax The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "GenLambdaMax", function(object) { c(1,1) })

#' @describeIn GenLambdaMax The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "GenLambdaMax", function(object) { c(FALSE, FALSE) })

#' @describeIn GenLambdaMax Is the atom convex?
setMethod("is_atom_convex", "GenLambdaMax", function(object) { FALSE })

#' @describeIn GenLambdaMax Is the atom concave?
setMethod("is_atom_concave", "GenLambdaMax", function(object) { FALSE })

#' @describeIn GenLambdaMax Is the atom quasiconvex?
setMethod("is_atom_quasiconvex", "GenLambdaMax", function(object) { TRUE })

#' @describeIn GenLambdaMax Is the atom quasiconcave?
setMethod("is_atom_quasiconcave", "GenLambdaMax", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn GenLambdaMax Is the atom weakly increasing in the index?
setMethod("is_incr", "GenLambdaMax", function(object, idx) { FALSE })

#' @describeIn GenLambdaMax Is the atom weakly decreasing in the index?
setMethod("is_decr", "GenLambdaMax", function(object, idx) { FALSE })

#'
#' The GeoMean class.
#'
#' This class represents the (weighted) geometric mean of vector \eqn{x} with optional powers given by \eqn{p}.
#'
#' \deqn{\left(x_1^{p_1} \cdots x_n^{p_n} \right)^{\frac{1}{\mathbf{1}^Tp}}}
#'
#' The geometric mean includes an implicit constraint that \eqn{x_i \geq 0} whenever \eqn{p_i > 0}. If \eqn{p_i = 0, x_i} will be unconstrained.
#' The only exception to this rule occurs when \eqn{p} has exactly one nonzero element, say \eqn{p_i}, in which case \code{GeoMean(x,p)} is equivalent to \eqn{x_i} (without the nonnegativity constraint).
#' A specific case of this is when \eqn{x \in \mathbf{R}^1}.
#'
#' @slot x An \linkS4class{Expression} or numeric vector.
#' @slot p (Optional) A vector of weights for the weighted geometric mean. The default is a vector of ones, giving the \strong{unweighted} geometric mean \eqn{x_1^{1/n} \cdots x_n^{1/n}}.
#' @slot max_denom (Optional) The maximum denominator to use in approximating \code{p/sum(p)} with \code{w}. If \code{w} is not an exact representation, increasing \code{max_denom} may offer a more accurate representation, at the cost of requiring more convex inequalities to represent the geometric mean. Defaults to 1024.
#' @slot w (Internal) A list of \code{bigq} objects that represent a rational approximation of \code{p/sum(p)}.
#' @slot approx_error (Internal) The error in approximating \code{p/sum(p)} with \code{w}, given by \eqn{\|p/\mathbf{1}^Tp - w\|_{\infty}}.
#' @name GeoMean-class
#' @aliases GeoMean
#' @importClassesFrom gmp bigq bigz
#' @rdname GeoMean-class
.GeoMean <- setClass("GeoMean", representation(x = "ConstValORExpr", p = "numeric", max_denom = "numeric",
                                               w = "bigq", w_dyad = "bigq", approx_error = "numeric", tree = "Rdict",
                                               cone_lb = "numeric", cone_num = "numeric", cone_num_over = "numeric"),
                                prototype(p = NA_real_, max_denom = 1024), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric vector.
#' @param p (Optional) A vector of weights for the weighted geometric mean. The default is a vector of ones, giving the \strong{unweighted} geometric mean \eqn{x_1^{1/n} \cdots x_n^{1/n}}.
#' @param max_denom (Optional) The maximum denominator to use in approximating \code{p/sum(p)} with \code{w}. If \code{w} is not an exact representation, increasing \code{max_denom} may offer a more accurate representation, at the cost of requiring more convex inequalities to represent the geometric mean. Defaults to 1024.
#' @rdname GeoMean-class
GeoMean <- function(x, p = NA_real_, max_denom = 1024) { .GeoMean(x = x, p = p, max_denom  = max_denom) }

setMethod("initialize", "GeoMean", function(.Object, ..., x, p = NA_real_, max_denom = 1024) {
  .Object@x <- x
  .Object@max_denom <- max_denom
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@x), validate = FALSE)

  x <- .Object@args[[1]]
  if(is_vector(x))
    n <- ifelse(ndim(x) == 0, 1, max(dim(x)))
  else
    stop("x must be a row or column vector.")

  if(any(is.na(p)))
    p <- rep(1, n)
  .Object@p <- p

  if(length(.Object@p) != n)
    stop("x and p must have the same number of elements.")

  if(any(.Object@p < 0) || sum(.Object@p) <= 0)
    stop("powers must be nonnegative and not all zero.")

  frac <- fracify(.Object@p, .Object@max_denom)
  .Object@w <- frac[[1]]
  .Object@w_dyad <- frac[[2]]
  .Object@approx_error <- approx_error(.Object@p, .Object@w)

  .Object@tree <- decompose(.Object@w_dyad)

  # known lower bound on number of cones needed to represent w_dyad
  .Object@cone_lb <- lower_bound(.Object@w_dyad)

  # number of cones used past known lower bound
  .Object@cone_num_over <- over_bound(.Object@w_dyad, .Object@tree)

  # number of cones used
  .Object@cone_num <- .Object@cone_lb + .Object@cone_num_over
  validObject(.Object)
  .Object
})

#' @param object A \linkS4class{GeoMean} object.
#' @param values A list of arguments to the atom.
#' @describeIn GeoMean The (weighted) geometric mean of the elements of \code{x}.
setMethod("to_numeric", "GeoMean", function(object, values) {
  values <- as.vector(values[[1]])
  val <- 1.0
  for(idx in 1:length(values)) {
    x <- values[[idx]]
    p <- object@w[idx]
    val <- val * Rmpfr::mpfr(x, Rmpfr::getPrec(x))^p
  }
  return(gmp::asNumeric(val))   # TODO: Handle mpfr objects in the backend later
})

#' @describeIn GeoMean Returns constraints describing the domain of the node
setMethod(".domain", "GeoMean", function(object) { list(object@args[[1]][object@w > 0] >= 0) })

#' @param values A list of numeric values for the arguments
#' @describeIn GeoMean Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "GeoMean", function(object, values) {
  x <- as.matrix(values[[1]])
  # No special case when only one non-zero weight
  w_arr <- as.double(object@w)   # TODO: I'm casting bigq/bigz to double to construct Matrix properly.
  # Outside domain
  if(any(x[w_arr > 0] <= 0))
    return(list(NA_real_))
  else {
    D <- w_arr/as.vector(x) * to_numeric(object, values)
    return(list(Matrix(D, sparse = TRUE)))
  }
})

#' @describeIn GeoMean The name and arguments of the atom.
setMethod("name", "GeoMean", function(x) {
  vals <- paste(sapply(x@w, as.character), collapse = ", ")
  paste("GeoMean(", name(x@args[[1]]), ", (", vals, "))", sep = "")
})

#' @describeIn GeoMean The atom is a scalar.
setMethod("dim_from_args", "GeoMean", function(object) { c(1,1) })

#' @describeIn GeoMean The atom is non-negative.
setMethod("sign_from_args", "GeoMean", function(object) { c(TRUE, FALSE) })

#' @describeIn GeoMean The atom is not convex.
setMethod("is_atom_convex", "GeoMean", function(object) { FALSE })

#' @describeIn GeoMean The atom is concave.
setMethod("is_atom_concave", "GeoMean", function(object) { TRUE })

#' @describeIn GeoMean Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "GeoMean", function(object) { TRUE })

#' @describeIn GeoMean Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "GeoMean", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn GeoMean The atom is weakly increasing in every argument.
setMethod("is_incr", "GeoMean", function(object, idx) { TRUE })

#' @describeIn GeoMean The atom is not weakly decreasing in any argument.
setMethod("is_decr", "GeoMean", function(object, idx) { FALSE })

#' @describeIn GeoMean Returns \code{list(w, dyadic completion, tree of dyads)}.
setMethod("get_data", "GeoMean", function(object) { list(object@w, object@w_dyad, object@tree) })

#' @param args An optional list that contains the arguments to reconstruct the atom. Default is to use current arguments of the atom.
#' @param id_objects Currently unused.
#' @describeIn GeoMean Returns a shallow copy of the GeoMean atom
setMethod("copy", "GeoMean", function(object, args = NULL, id_objects = list()) {
  if(is.null(args))
    args <- object@args

  copy <- do.call(class(object), args)
  data <- get_data(object)
  copy@w <- data[[1]]
  copy@w_dyad <- data[[2]]
  copy@tree <- data[[3]]

  copy@approx_error <- object@approx_error
  copy@cone_lb <- object@cone_lb
  copy@cone_num_over <- object@cone_num_over
  copy@cone_num <- object@cone_num
  copy
})

#'
#' The HarmonicMean atom.
#'
#' The harmonic mean of \eqn{x}, \eqn{\frac{1}{n} \sum_{i=1}^n x_i^{-1}}, where \eqn{n} is the length of \eqn{x}.
#'
#' @param x An expression or number whose harmonic mean is to be computed. Must have positive entries.
#' @return The harmonic mean of \code{x}.
HarmonicMean <- function(x) {
  x <- as.Constant(x)
  size(x) * Pnorm(x = x, p = -1)
}

#'
#' The InvProb atom.
#' 
#' The reciprocal of a product of the entries of a vector \eqn{x}, \eqn{(\prod_{i=1}^n x_i)^{-1}}, where \eqn{n} is the length of \eqn{x}.
#'
#' @param x An expression or vector whose reciprocal product is to be computed. Must have positive entries.
#' @return The reciprocal product of \code{x}.
InvProb <- function(value) {
  val_dim <- dim(value)
  if(is.na(val_dim) || is.null(val_dim))
    p <- 1
  else
    p <- as.integer(sum(val_dim))
  Power(InvPos(GeoMean(value)), p)
}

#'
#' The LambdaMax class.
#'
#' The maximum eigenvalue of a matrix, \eqn{\lambda_{\max}(A)}.
#'
#' @slot A An \linkS4class{Expression} or numeric matrix.
#' @name LambdaMax-class
#' @aliases LambdaMax
#' @rdname LambdaMax-class
.LambdaMax <- setClass("LambdaMax", representation(A = "ConstValORExpr"), contains = "Atom")

#' @param A An \linkS4class{Expression} or numeric matrix.
#' @rdname LambdaMax-class
LambdaMax <- function(A) { .LambdaMax(A = A) }

setMethod("initialize", "LambdaMax", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., atom_args = list(.Object@A))
})

#' @param object A \linkS4class{LambdaMax} object.
#' @param values A list of arguments to the atom.
#' @describeIn LambdaMax The largest eigenvalue of \code{A}. Requires that \code{A} be symmetric.
setMethod("to_numeric", "LambdaMax", function(object, values) {
  # if(any(t(values[[1]]) != values[[1]]))
  #  stop("LambdaMax called on a non-symmetric matrix")
  max(eigen(values[[1]], only.values = TRUE)$values)
})

#' @describeIn LambdaMax Returns the constraints describing the domain of the atom.
setMethod(".domain", "LambdaMax", function(object) { list(Conj(t(object@args[[1]])) == object@args[[1]]) })

#' @describeIn LambdaMax Gives the (sub/super)gradient of the atom with respect to each argument. Matrix expressions are vectorized, so the gradient is a matrix.
setMethod(".grad", "LambdaMax", function(object, values) {
  r <- base::eigen(values[[1]], only.values = FALSE)   # Eigenvalues returned in decreasing order.
  v <- r$vectors  # eigenvectors
  w <- r$values   # eigenvalues

  d <- rep(0, length(w))
  d[1] <- 1
  d <- diag(d)
  D <- v %*% d %*% t(v)
  list(Matrix(as.vector(D), sparse = TRUE))
})

#' @describeIn LambdaMax Check that \code{A} is square.
setMethod("validate_args", "LambdaMax", function(object) {
  if(ndim(object@args[[1]]) != 2 || nrow(object@args[[1]]) != ncol(object@args[[1]]))
    stop("The argument to LambdaMax must resolve to a square matrix")
})

#' @describeIn LambdaMax The atom is a scalar.
setMethod("dim_from_args", "LambdaMax", function(object) { c(1,1) })

#' @describeIn LambdaMax The sign of the atom is unknown.
setMethod("sign_from_args", "LambdaMax", function(object) { c(FALSE, FALSE) })

#' @describeIn LambdaMax The atom is convex.
setMethod("is_atom_convex", "LambdaMax", function(object) { TRUE })

#' @describeIn LambdaMax The atom is not concave.
setMethod("is_atom_concave", "LambdaMax", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn LambdaMax The atom is not monotonic in any argument.
setMethod("is_incr", "LambdaMax", function(object, idx) { FALSE })

#' @describeIn LambdaMax The atom is not monotonic in any argument.
setMethod("is_decr", "LambdaMax", function(object, idx) { FALSE })

#'
#' The LambdaMin atom.
#'
#' The minimum eigenvalue of a matrix, \eqn{\lambda_{\min}(A)}.
#'
#' @param A An \linkS4class{Expression} or numeric matrix.
#' @return Returns the minimum eigenvalue of a matrix.
LambdaMin <- function(A) {
  A <- as.Constant(A)
  -LambdaMax(-A)
}

#'
#' The LambdaSumLargest class.
#'
#' This class represents the sum of the \code{k} largest eigenvalues of a matrix.
#'
#' @slot k A positive integer.
#' @name LambdaSumLargest-class
#' @aliases LambdaSumLargest
#' @rdname LambdaSumLargest-class
.LambdaSumLargest <- setClass("LambdaSumLargest", representation(k = "numeric"), contains = "LambdaMax")

#' @param A An \linkS4class{Expression} or numeric matrix.
#' @param k A positive integer.
#' @rdname LambdaSumLargest-class
LambdaSumLargest <- function(A, k) { .LambdaSumLargest(A = A, k = k) }

setMethod("initialize", "LambdaSumLargest", function(.Object, ..., k) {
  .Object@k <- k
  callNextMethod(.Object, ...)
})

#' @describeIn LambdaSumLargest Does the atom handle complex numbers?
setMethod("allow_complex", "LambdaSumLargest", function(object) { TRUE })

#' @param object A \linkS4class{LambdaSumLargest} object.
#' @param values A list of arguments to the atom.
#' @describeIn LambdaSumLargest Returns the largest eigenvalue of \code{A}, which must be symmetric.
setMethod("to_numeric", "LambdaSumLargest", function(object, values) {
  # if(any(t(values[[1]]) != values[[1]]))
  #  stop("LambdaSumLargest called on a non-symmetric matrix")
  eigs <- eigen(values[[1]], only.values = TRUE)$values
  value(SumLargest(eigs, object@k))
})

#' @describeIn LambdaSumLargest Verify that the argument \code{A} is square.
setMethod("validate_args", "LambdaSumLargest", function(object) {
  A <- object@args[[1]]
  if(ndim(A) != 2 || nrow(A) != ncol(A))
    stop("First argument must be a square matrix.")
  else if(as.integer(object@k) != object@k || object@k <= 0)
    stop("Second argument must be a positive integer.")
})

#' @describeIn LambdaSumLargest Returns the parameter \code{k}.
setMethod("get_data", "LambdaSumLargest", function(object) { list(object@k) })

#' @param values A list of numeric values for the arguments
#' @describeIn LambdaSumLargest Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "LambdaSumLargest", function(object, values) { stop("Unimplemented") })

#'
#' The LambdaSumSmallest atom.
#'
#' This class represents the sum of the \code{k} smallest eigenvalues of a matrix.
#'
#' @param A An \linkS4class{Expression} or numeric matrix.
#' @param k A positive integer.
#' @return Returns the sum of the k smallest eigenvalues of a matrix.
LambdaSumSmallest <- function(A, k) {
  A <- as.Constant(A)
  -LambdaSumLargest(-A, k)
}

#'
#' The VecLength class.
#'
#' This class represents the length of a vector (index of last nonzero, ones-based).
#'
#' @slot x An \linkS4class{Expression} representing a vector.
#' @name VecLength-class
#' @aliases VecLength
#' @rdname VecLength-class
.VecLength <- setClass("VecLength", representation(x = "ConstValORExpr"),
                         validity = function(object) {
                           if(!is_vector(object@args[[1]]))
                             stop("[VecLength: x] The argument x must be a vector.")
                           return(TRUE)
                         }, contains = "Atom")

#' @param x An \linkS4class{Expression} representing a vector.
#' @rdname VecLength-class
VecLength <- function(x) { .Length(x = x) }

setMethod("initialize", "VecLength", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object A \linkS4class{VecLength} object.
#' @param values A list of arguments to the atom.
#' @describeIn VecLength The length of the vector.
setMethod("to_numeric", "VecLength", function(object, values) {
  outside_tol <- abs(values[[1]]) > ATOM_EVAL_TOL
  return(base::max(which(outside_tol)))
})

#' @describeIn VecLength The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "VecLength", function(object) { c(1,1) })

#' @describeIn VecLength The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "VecLength", function(object) { c(TRUE, FALSE) })

#' @describeIn VecLength Is the atom convex?
setMethod("is_atom_convex", "VecLength", function(object) { FALSE })

#' @describeIn VecLength Is the atom concave?
setMethod("is_atom_concave", "VecLength", function(object) { FALSE })

#' @describeIn VecLength Is the atom quasiconvex?
setMethod("is_atom_quasiconvex", "VecLength", function(object) { TRUE })

#' @describeIn VecLength Is the atom quasiconcave?
setMethod("is_atom_quasiconcave", "VecLength", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn VecLength Is the atom weakly increasing in the index?
setMethod("is_incr", "VecLength", function(object, idx) { FALSE })

#' @describeIn VecLength Is the atom weakly decreasing in the index?
setMethod("is_decr", "VecLength", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn VecLength Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "VecLength", function(object, values) { NA_real_ })

#'
#' The LogDet class.
#'
#' The natural logarithm of the determinant of a matrix, \eqn{\log\det(A)}.
#'
#' @slot A An \linkS4class{Expression} or numeric matrix.
#' @name LogDet-class
#' @aliases LogDet
#' @rdname LogDet-class
.LogDet <- setClass("LogDet", representation(A = "ConstValORExpr"), contains = "Atom")

#' @param A An \linkS4class{Expression} or numeric matrix.
#' @rdname LogDet-class
LogDet <- function(A) { .LogDet(A = A) }

setMethod("initialize", "LogDet", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., atom_args = list(.Object@A))
})

#' @param object A \linkS4class{LogDet} object.
#' @param values A list of arguments to the atom.
#' @describeIn LogDet The log-determinant of SDP matrix \code{A}. This is the sum of logs of the eigenvalues and is equivalent to the nuclear norm of the matrix logarithm of \code{A}.
setMethod("to_numeric", "LogDet", function(object, values) {
  if(is.complex(values[[1]])) {
    eigvals <- eigen(values[[1]], only.values = TRUE)$values
    return(log(prod(eigvals)))
  } else {
    logdet <- determinant(values[[1]], logarithm = TRUE)
    if(logdet$sign == 1)
      return(as.numeric(logdet$modulus))
    else
      return(-Inf)
  }
})

#' @describeIn LogDet Check that \code{A} is square.
setMethod("validate_args", "LogDet", function(object) {
  arg_dim <- dim(object@args[[1]])
  if(length(arg_dim) == 1 || arg_dim[1] != arg_dim[2])
    stop("The argument to LogDet must be a square matrix")
})

#' @describeIn LogDet The atom is a scalar.
setMethod("dim_from_args", "LogDet", function(object) { c(1,1) })

#' @describeIn LogDet The atom is non-negative.
setMethod("sign_from_args",  "LogDet", function(object) { c(TRUE, FALSE) })

#' @describeIn LogDet The atom is not convex.
setMethod("is_atom_convex", "LogDet", function(object) { FALSE })

#' @describeIn LogDet The atom is concave.
setMethod("is_atom_concave", "LogDet", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn LogDet The atom is not monotonic in any argument.
setMethod("is_incr", "LogDet", function(object, idx) { FALSE })

#' @describeIn LogDet The atom is not monotonic in any argument.
setMethod("is_decr", "LogDet", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn LogDet Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "LogDet", function(object, values) {
  X <- as.matrix(values[[1]])
  eigen_val <- eigen(X, only.values = TRUE)$values
  if(min(eigen_val) > 0) {
    # Grad: t(X^(-1))
    D <- t(base::solve(X))
    return(list(Matrix(as.vector(D), sparse = TRUE)))
  } else   # Outside domain
    return(list(NA_real_))
})

#' @describeIn LogDet Returns constraints describing the domain of the node
setMethod(".domain", "LogDet", function(object) { list(object@args[[1]] %>>% 0) })

#'
#' The LogSumExp class.
#'
#' The natural logarithm of the sum of the elementwise exponential, \eqn{\log\sum_{i=1}^n e^{x_i}}.
#'
#' @slot x An \linkS4class{Expression} representing a vector or matrix.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @slot keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @name LogSumExp-class
#' @aliases LogSumExp
#' @rdname LogSumExp-class
.LogSumExp <- setClass("LogSumExp", contains = "AxisAtom")

#' @param x An \linkS4class{Expression} representing a vector or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @rdname LogSumExp-class
LogSumExp <- function(x, axis = NA_real_, keepdims = FALSE) { .LogSumExp(expr = x, axis = axis, keepdims = keepdims) }

#' @param object A \linkS4class{LogSumExp} object.
#' @param values A list of arguments to the atom.
#' @describeIn LogSumExp Evaluates \eqn{e^x} elementwise, sums, and takes the natural log.
setMethod("to_numeric", "LogSumExp", function(object, values) {
  if(is.na(object@axis))
    log(sum(exp(values[[1]])))
  else
    # log(apply(exp(values[[1]]), object@axis, sum))
    log(apply_with_keepdims(exp(values[[1]]), sum, axis = object@axis, keepdims = object@keepdims))
})

#' @param values A list of numeric values.
#' @describeIn LogSumExp Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "LogSumExp", function(object, values) { .axis_grad(object, values) })

#' @param value A numeric value.
#' @describeIn LogSumExp Gives the (sub/super)gradient of the atom w.r.t. each column variable.
setMethod(".column_grad", "LogSumExp", function(object, value) {
  denom <- sum(exp(value))
  nom <- exp(value)
  D <- nom/denom
  D
})

#' @describeIn LogSumExp Returns sign (is positive, is negative) of the atom.
setMethod("sign_from_args",  "LogSumExp", function(object) { c(FALSE, FALSE) })

#' @describeIn LogSumExp The atom is convex.
setMethod("is_atom_convex", "LogSumExp", function(object) { TRUE })

#' @describeIn LogSumExp The atom is not concave.
setMethod("is_atom_concave", "LogSumExp", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn LogSumExp The atom is weakly increasing in the index.
setMethod("is_incr", "LogSumExp", function(object, idx) { TRUE })

#' @param idx An index into the atom.
#' @describeIn LogSumExp The atom is not weakly decreasing in the index.
setMethod("is_decr", "LogSumExp", function(object, idx) { FALSE })

#'
#' The MatrixFrac class.
#'
#' The matrix fraction function \eqn{tr(X^T P^{-1} X)}.
#'
#' @slot X An \linkS4class{Expression} or numeric matrix.
#' @slot P An \linkS4class{Expression} or numeric matrix.
#' @name MatrixFrac-class
#' @aliases MatrixFrac
#' @rdname MatrixFrac-class
.MatrixFrac <- setClass("MatrixFrac", representation(X = "ConstValORExpr", P = "ConstValORExpr"), contains = "Atom")

#' @param X An \linkS4class{Expression} or numeric matrix.
#' @param P An \linkS4class{Expression} or numeric matrix.
#' @rdname MatrixFrac-class
MatrixFrac <- function(X, P) { .MatrixFrac(X = X, P = P) }

setMethod("initialize", "MatrixFrac", function(.Object, ..., X, P) {
  .Object@X <- X
  .Object@P <- P
  callNextMethod(.Object, ..., atom_args = list(.Object@X, .Object@P))
})

#' @describeIn MatrixFrac Does the atom handle complex numbers?
setMethod("allow_complex", "MatrixFrac", function(object) { TRUE })

#' @param object A \linkS4class{MatrixFrac} object.
#' @param values A list of arguments to the atom.
#' @describeIn MatrixFrac The trace of \eqn{X^TP^{-1}X}.
setMethod("to_numeric", "MatrixFrac", function(object, values) {
  # TODO: Raise error if not invertible?
  X <- values[[1]]
  P <- values[[2]]
  if(is_complex(object@args[[1]]))
    product <- t(Conj(X)) %*% base::solve(P) %*% X
  else
    product <- t(X) %*% base::solve(P) %*% X

  if(length(dim(product)) == 2)
    return(sum(diag(product)))
  else
    return(product)
})

#' @describeIn MatrixFrac Check that the dimensions of \code{x} and \code{P} match.
setMethod("validate_args", "MatrixFrac", function(object) {
  X <- object@args[[1]]
  P <- object@args[[2]]
  if(ndim(P) != 2 || nrow(P) != ncol(P))
    stop("The second argument to MatrixFrac must be a square matrix.")
  else if(nrow(X) != nrow(P))
    stop("The arguments to MatrixFrac have incompatible dimensions.")
})

#' @describeIn MatrixFrac The atom is a scalar.
setMethod("dim_from_args", "MatrixFrac", function(object) { c(1,1) })

#' @describeIn MatrixFrac The atom is positive.
setMethod("sign_from_args", "MatrixFrac", function(object) { c(TRUE, FALSE) })

#' @describeIn MatrixFrac The atom is convex.
setMethod("is_atom_convex", "MatrixFrac", function(object) { TRUE })

#' @describeIn MatrixFrac The atom is not concave.
setMethod("is_atom_concave", "MatrixFrac", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn MatrixFrac The atom is not monotonic in any argument.
setMethod("is_incr", "MatrixFrac", function(object, idx) { FALSE })

#' @describeIn MatrixFrac The atom is not monotonic in any argument.
setMethod("is_decr", "MatrixFrac", function(object, idx) { FALSE })

#' @describeIn MatrixFrac True if x is affine and P is constant.
setMethod("is_quadratic", "MatrixFrac", function(object) { is_affine(object@args[[1]]) && is_constant(object@args[[2]]) })

#' @describeIn MatrixFrac True if x is piecewise linear and P is constant.
setMethod("is_qpwa", "MatrixFrac", function(object) { is_pwl(object@args[[1]]) && is_constant(object@args[[2]]) })

#' @describeIn MatrixFrac Returns constraints describing the domain of the node
setMethod(".domain", "MatrixFrac", function(object) { list(object@args[[2]] %>>% 0) })

#' @param values A list of numeric values for the arguments
#' @describeIn MatrixFrac Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "MatrixFrac", function(object, values) {
  X <- as.matrix(values[[1]])
  P <- as.matrix(values[[2]])
  P_inv <- tryCatch({
    base::solve(P)
  }, error = function(e) {
      list(NA_real_, NA_real_)
  })

  ## if(is.null(dim(P_inv)) && is.na(P_inv))
  ##   return(list(NA_real_, NA_real_))

  if(is.null(dim(P_inv)))
      return(list(NA_real_, NA_real_))

  # partial_X = (P^-1+P^-T)X
  # partial_P = (P^-1 * X * X^T * P^-1)^T
  DX <- (P_inv + t(P_inv)) %*% X
  DX <- as.vector(t(DX))
  DX <- Matrix(DX, sparse = TRUE)

  DP <- P_inv %*% X
  DP <- DP %*% t(X)
  DP <- DP %*% P_inv
  DP <- -t(DP)
  DP <- Matrix(as.vector(t(DP)), sparse = TRUE)
  list(DX, DP)
})

#'
#' The MaxEntries class.
#'
#' The maximum of an expression.
#'
#' @slot x An \linkS4class{Expression} representing a vector or matrix.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @slot keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @name MaxEntries-class
#' @aliases MaxEntries
#' @rdname MaxEntries-class
.MaxEntries <- setClass("MaxEntries", contains = "AxisAtom")

#' @param x An \linkS4class{Expression} representing a vector or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @rdname MaxEntries-class
MaxEntries <- function(x, axis = NA_real_, keepdims = FALSE) { .MaxEntries(expr = x, axis = axis, keepdims = keepdims) }

#' @param object A \linkS4class{MaxEntries} object.
#' @param values A list of arguments to the atom.
#' @describeIn MaxEntries The largest entry in \code{x}.
setMethod("to_numeric", "MaxEntries", function(object, values) {
  apply_with_keepdims(values[[1]], max, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn MaxEntries The sign of the atom.
setMethod("sign_from_args",  "MaxEntries", function(object) { c(is_nonneg(object@args[[1]]), is_nonpos(object@args[[1]])) })

#' @describeIn MaxEntries The atom is convex.
setMethod("is_atom_convex", "MaxEntries", function(object) { TRUE })

#' @describeIn MaxEntries The atom is not concave.
setMethod("is_atom_concave", "MaxEntries", function(object) { FALSE })

#' @describeIn MaxEntries Is the atom log-log convex.
setMethod("is_atom_log_log_convex", "MaxEntries", function(object) { TRUE })

#' @describeIn MaxEntries Is the atom log-log concave.
setMethod("is_atom_log_log_concave", "MaxEntries", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn MaxEntries The atom is weakly increasing in every argument.
setMethod("is_incr", "MaxEntries", function(object, idx) { TRUE })

#' @param idx An index into the atom.
#' @describeIn MaxEntries The atom is not weakly decreasing in any argument.
setMethod("is_decr", "MaxEntries", function(object, idx) { FALSE })

#' @describeIn MaxEntries Is \code{x} piecewise linear?
setMethod("is_pwl", "MaxEntries", function(object) { is_pwl(object@args[[1]]) })

#' @param values A list of numeric values for the arguments
#' @describeIn MaxEntries Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "MaxEntries", function(object, values) { .axis_grad(object, values) })

#' @param value A numeric value
#' @describeIn MaxEntries Gives the (sub/super)gradient of the atom w.r.t. each column variable
setMethod(".column_grad", "MaxEntries", function(object, value) {
  # Grad: 1 for a largest index
  value <- as.vector(value)
  idx <- (value == max(value))
  D <- matrix(0, nrow = length(value), ncol = 1)
  D[idx,1] <- 1
  D
})

#'
#' The MinEntries class.
#'
#' The minimum of an expression.
#'
#' @slot x An \linkS4class{Expression} representing a vector or matrix.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @slot keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @name MinEntries-class
#' @aliases MinEntries
#' @rdname MinEntries-class
.MinEntries <- setClass("MinEntries", contains = "AxisAtom")

#' @param x An \linkS4class{Expression} representing a vector or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @rdname MinEntries-class
MinEntries <- function(x, axis = NA_real_, keepdims = FALSE) { .MinEntries(expr = x, axis = axis, keepdims = keepdims) }

#' @param object A \linkS4class{MinEntries} object.
#' @param values A list of arguments to the atom.
#' @describeIn MinEntries The largest entry in \code{x}.
setMethod("to_numeric", "MinEntries", function(object, values) {
  apply_with_keepdims(values[[1]], min, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn MinEntries The sign of the atom.
setMethod("sign_from_args",  "MinEntries", function(object) { c(is_nonneg(object@args[[1]]), is_nonpos(object@args[[1]])) })

#' @describeIn MinEntries The atom is not convex.
setMethod("is_atom_convex", "MinEntries", function(object) { FALSE })

#' @describeIn MinEntries The atom is concave.
setMethod("is_atom_concave", "MinEntries", function(object) { TRUE })

#' @describeIn MinEntries Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "MinEntries", function(object) { FALSE })

#' @describeIn MinEntries Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "MinEntries", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn MinEntries The atom is weakly increasing in every argument.
setMethod("is_incr", "MinEntries", function(object, idx) { TRUE })

#' @param idx An index into the atom.
#' @describeIn MinEntries The atom is not weakly decreasing in any argument.
setMethod("is_decr", "MinEntries", function(object, idx) { FALSE })

#' @describeIn MinEntries Is \code{x} piecewise linear?
setMethod("is_pwl", "MinEntries", function(object) { is_pwl(object@args[[1]]) })

#' @param values A list of numeric values for the arguments
#' @describeIn MinEntries Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "MinEntries", function(object, values) { .axis_grad(object, values) })

#' @param value A numeric value
#' @describeIn MinEntries Gives the (sub/super)gradient of the atom w.r.t. each column variable
setMethod(".column_grad", "MinEntries", function(object, value) {
  # Grad: 1 for a largest index
  value <- as.vector(value)
  idx <- (value == min(value))
  D <- matrix(0, nrow = length(value), ncol = 1)
  D[idx,1] <- 1
  D
})

#'
#' The Pnorm class.
#'
#' This class represents the vector p-norm.
#'
#' If given a matrix variable, \code{Pnorm} will treat it as a vector and compute the p-norm of the concatenated columns.
#'
#' For \eqn{p \geq 1}, the p-norm is given by \deqn{\|x\|_p = \left(\sum_{i=1}^n |x_i|^p\right)^{1/p}} with domain \eqn{x \in \mathbf{R}^n}.
#' For \eqn{p < 1, p\neq 0}, the p-norm is given by \deqn{\|x\|_p = \left(\sum_{i=1}^n x_i^p\right)^{1/p}} with domain \eqn{x \in \mathbf{R}^n_+}.
#'
#' \itemize{
#'    \item Note that the "p-norm" is actually a \strong{norm} only when \eqn{p \geq 1} or \eqn{p = +\infty}. For these cases, it is convex.
#'    \item The expression is undefined when \eqn{p = 0}.
#'    \item Otherwise, when \eqn{p < 1}, the expression is concave, but not a true norm.
#' }
#'
#' @slot x An \linkS4class{Expression} representing a vector or matrix.
#' @slot p A number greater than or equal to 1, or equal to positive infinity.
#' @slot max_denom The maximum denominator considered in forming a rational approximation for \eqn{p}.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @slot keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @slot .approx_error (Internal) The absolute difference between \eqn{p} and its rational approximation.
#' @slot .original_p (Internal) The original input \eqn{p}.
#' @name Pnorm-class
#' @aliases Pnorm
#' @rdname Pnorm-class
.Pnorm <- setClass("Pnorm", representation(p = "numeric", max_denom = "numeric", .approx_error = "numeric", .original_p = "numeric"),
                  prototype(p = 2, max_denom = 1024, .approx_error = NA_real_, .original_p = NA_real_), contains = "AxisAtom")

#' @param x An \linkS4class{Expression} representing a vector or matrix.
#' @param p A number greater than or equal to 1, or equal to positive infinity.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @param max_denom (Optional) The maximum denominator considered in forming a rational approximation for \eqn{p}. The default is 1024.
#' @rdname Pnorm-class
Pnorm <- function(x, p = 2, axis = NA_real_, keepdims = FALSE, max_denom = 1024) {
  if(p == 1)
    Norm1(x, axis = axis, keepdims = keepdims)
  else if(p %in% c(Inf, "inf", "Inf"))
    NormInf(x, axis = axis, keepdims = keepdims)
  else
    .Pnorm(expr = x, axis = axis, keepdims = keepdims, p = p, max_denom = max_denom)
}

setMethod("initialize", "Pnorm", function(.Object, ..., p = 2, max_denom = 1024, .approx_error = NA_real_, .original_p = NA_real_) {
  if(p == 1)
    stop("Use the Norm1 class to instantiate a 1-norm.")
  else if(p %in% c(Inf, "inf", "Inf"))
    stop("Use the NormInf class to instantiate an infinity-norm.")
  # else if(p < 0)
  #  .Object@p <- pow_neg(p, max_denom)
  # else if(p > 0 && p < 1)
  #  .Object@p <- pow_mid(p, max_denom)
  # else if(p > 1)
  #  .Object@p <- pow_high(p, max_denom)
  # else
  #  stop("Invalid value of p.")
  .Object@p <- p

  .Object@max_denom <- max_denom
  .Object@.approx_error <- abs(.Object@p - p)
  .Object@.original_p <- p
  callNextMethod(.Object, ...)
})

#'
#' Internal method for calculating the p-norm
#'
#' @param x A matrix
#' @param p A number grater than or equal to 1, or equal to positive infinity
#' @return Returns the specified norm of matrix x
.p_norm <- function(x, p) {
  if(p == Inf)
    max(abs(x))
  else if(p == 0)
    sum(x != 0)
  else if(p %% 2 == 0 || p < 1)
    sum(x^p)^(1/p)
  else if(p >= 1)
    sum(abs(x)^p)^(1/p)
  else
    stop("Invalid p = ", p)
}

#' @describeIn Pnorm Does the atom handle complex numbers?
setMethod("allow_complex", "Pnorm", function(object) { TRUE })

#' @param object A \linkS4class{Pnorm} object.
#' @param values A list of arguments to the atom.
#' @describeIn Pnorm The p-norm of \code{x}.
setMethod("to_numeric", "Pnorm", function(object, values) {
  if(is.na(object@axis))
    values <- as.vector(values[[1]])
  else
    values <- as.matrix(values[[1]])

  if(object@p < 1 && any(values < 0))
    return(-Inf)
  if(object@p < 0 && any(values == 0))
    return(0)

  apply_with_keepdims(values, function(x) { .p_norm(x, object@p) }, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn Pnorm Check that the arguments are valid.
setMethod("validate_args", "Pnorm", function(object) {
  callNextMethod()
  if(!is.na(object@axis) && object@p != 2)
    stop("The axis parameter is only supported for p = 2.")
  if(object@p < 1 && is_complex(object@args[[1]]))
    stop("Pnorm(x, p) cannot have x complex for p < 1.")
})

#' @describeIn Pnorm The atom is positive.
setMethod("sign_from_args",  "Pnorm", function(object) { c(TRUE, FALSE) })

#' @describeIn Pnorm The atom is convex if \eqn{p \geq 1}.
setMethod("is_atom_convex", "Pnorm", function(object) { object@p > 1 })

#' @describeIn Pnorm The atom is concave if \eqn{p < 1}.
setMethod("is_atom_concave", "Pnorm", function(object) { object@p < 1 })

#' @describeIn Pnorm Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Pnorm", function(object) { TRUE })

#' @describeIn Pnorm Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Pnorm", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn Pnorm The atom is weakly increasing if \eqn{p < 1} or \eqn{p > 1} and \code{x} is positive.
setMethod("is_incr", "Pnorm", function(object, idx) { object@p < 1 || (object@p > 1 && is_nonneg(object@args[[1]])) })

#' @param idx An index into the atom.
#' @describeIn Pnorm The atom is weakly decreasing if \eqn{p > 1} and \code{x} is negative.
setMethod("is_decr", "Pnorm", function(object, idx) { object@p > 1 && is_nonpos(object@args[[1]]) })

#' @describeIn Pnorm The atom is not piecewise linear unless \eqn{p = 1} or \eqn{p = \infty}.
setMethod("is_pwl", "Pnorm", function(object) { FALSE })

#' @describeIn Pnorm Returns \code{list(p, axis)}.
setMethod("get_data", "Pnorm", function(object) { list(object@p, object@axis) })

#' @describeIn Pnorm The name and arguments of the atom.
setMethod("name", "Pnorm", function(x) {
  sprintf("%s(%s, %s)", class(x), name(x@args[[1]]), x@p)
})

#' @describeIn Pnorm Returns constraints describing the domain of the node
setMethod(".domain", "Pnorm", function(object) {
  if(object@p < 1 && object@p != 0)
    list(object@args[[1]] >= 0)
  else
    list()
})

#' @param values A list of numeric values for the arguments
#' @describeIn Pnorm Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "Pnorm", function(object, values) { .axis_grad(object, values) })

#' @param value A numeric value
#' @describeIn Pnorm Gives the (sub/super)gradient of the atom w.r.t. each column variable
setMethod(".column_grad", "Pnorm", function(object, value) {
  rows <- size(object@args[[1]])
  value <- as.matrix(value)

  # Outside domain
  if(object@p < 1 && any(value <= 0))
    return(NA_real_)

  D_null <- sparseMatrix(i = c(), j = c(), dims = c(rows, 1))
  denominator <- .p_norm(value, object@p)
  denominator <- denominator^(object@p - 1)

  # Subgrad is 0 when denom is 0 (or undefined)
  if(denominator == 0) {
    if(object@p > 1)
      return(D_null)
    else
      return(NA_real_)
  } else {
    numerator <- value^(object@p - 1)
    frac <- numerator / denominator
    return(matrix(as.vector(frac)))
  }
})

#'
#' The MixedNorm atom.
#'
#' The \eqn{l_{p,q}} norm of X, \eqn{(\sum_k (\sum_l ||X_{k,l}||^p)^{q/p})^{1/q}}.
#'
#' @param X The matrix to take the \eqn{l_{p,q}} norm of
#' @param p The type of inner norm
#' @param q The type of outer norm
#' @return Returns the mixed norm of X with specified parameters p and q
MixedNorm <- function(X, p = 2, q = 1) {
  X <- as.Constant(X)

  # Inner norms
  vecnorms <- Norm(X, p, axis = 1)

  # Outer norms
  Norm(vecnorms, q)
}

#'
#' The Norm atom.
#'
#' Wrapper around the different norm atoms.
#'
#' @param x The matrix to take the norm of
#' @param p The type of norm. Valid options include any positive integer, 'fro' (for frobenius),
#' 'nuc' (sum of singular values), np.inf or 'inf' (infinity norm).
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @return Returns the specified norm of x.
#' @rdname Norm-atom
Norm <- function(x, p = 2, axis = NA_real_, keepdims = FALSE) {
  x <- as.Constant(x)

  # Matrix norms take precedence.
  num_nontrivial_idxs <- sum(dim(x) > 1)
  if(is.na(axis) && ndim(x) == 2) {
    if(p == 1)   # Matrix 1-norm.
      MaxEntries(Norm1(x, axis = 2))
    else if(p == "fro" || (p == 2 && num_nontrivial_idxs == 1))   # Frobenius norm.
      Pnorm(Vec(x), 2)
    else if(p == 2)   # Matrix 2-norm is largest singular value.
      SigmaMax(x)
    else if(p == "nuc")   # The nuclear norm (sum of singular values)
      NormNuc(x)
    else if(p %in% c(Inf, "inf", "Inf"))   # The matrix infinity-norm.
      MaxEntries(Norm1(x, axis = 1))
    else
      stop("Unsupported matrix norm.")
  } else {
    if(p == 1 || is_scalar(x))
      Norm1(x, axis = axis, keepdims = keepdims)
    else if(p %in% c(Inf, "inf", "Inf"))
      NormInf(x, axis = axis, keepdims = keepdims)
    else
      Pnorm(x, p, axis = axis, keepdims = keepdims)
  }
}

#'
#' The Norm2 atom.
#'
#' The 2-norm of an expression.
#'
#' @param x An \linkS4class{Expression} object.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @return Returns the 2-norm of x.
#' @rdname Norm2-atom
Norm2 <- function(x, axis = NA_real_, keepdims = FALSE) {
  Pnorm(x, p = 2, axis = axis, keepdims = keepdims)
}

#'
#' The Norm1 class.
#'
#' This class represents the 1-norm of an expression.
#'
#' @slot x An \linkS4class{Expression} object.
#' @name Norm1-class
#' @aliases Norm1
#' @rdname Norm1-class
.Norm1 <- setClass("Norm1", contains = "AxisAtom")

#' @param x An \linkS4class{Expression} object.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @rdname Norm1-class
Norm1 <- function(x, axis = NA_real_, keepdims = FALSE) { .Norm1(expr = x, axis = axis, keepdims = keepdims) }

#' @param object A \linkS4class{Norm1} object.
#' @describeIn Norm1 The name and arguments of the atom.
setMethod("name", "Norm1", function(x) {
  paste(class(x), "(", name(x@args[[1]]), ")", sep = "")
})

#' @param values A list of arguments to the atom.
#' @describeIn Norm1 Returns the 1-norm of x along the given axis.
setMethod("to_numeric", "Norm1", function(object, values) {
  if(is.na(object@axis))
    # base::norm(values[[1]], type = "O")
    sum(abs(values[[1]]))
  else
    # apply_with_keepdims(values[[1]], function(x) { norm(as.matrix(x), type = "O") }, axis = object@axis, keepdims = object@keepdims)
    apply_with_keepdims(values[[1]], function(x) { sum(abs(x)) }, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn Norm1 Does the atom handle complex numbers?
setMethod("allow_complex", "Norm1", function(object) { TRUE })

#' @describeIn Norm1 The atom is always positive.
setMethod("sign_from_args", "Norm1", function(object) { c(TRUE, FALSE) })

#' @describeIn Norm1 The atom is convex.
setMethod("is_atom_convex", "Norm1", function(object) { TRUE })

#' @describeIn Norm1 The atom is not concave.
setMethod("is_atom_concave", "Norm1", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn Norm1 Is the composition weakly increasing in argument \code{idx}?
setMethod("is_incr", "Norm1", function(object, idx) { is_nonneg(object@args[[1]]) })

#' @param idx An index into the atom.
#' @describeIn Norm1 Is the composition weakly decreasing in argument \code{idx}?
setMethod("is_decr", "Norm1", function(object, idx) { is_nonpos(object@args[[1]]) })

#' @describeIn Norm1 Is the atom piecewise linear?
setMethod("is_pwl", "Norm1", function(object) {
  is_pwl(object@args[[1]]) && (is_real(object@args[[1]]) || is_imag(object@args[[1]]))
})

#' @describeIn Norm1 Returns the axis.
setMethod("get_data", "Norm1", function(object) { list(object@axis) })

#' @describeIn Norm1 Returns constraints describing the domain of the node
setMethod(".domain", "Norm1", function(object) { list() })

#' @param values A list of numeric values for the arguments
#' @describeIn Norm1 Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "Norm1", function(object, values) { .axis_grad(object, values) })

#' @param value A numeric value
#' @describeIn Norm1 Gives the (sub/super)gradient of the atom w.r.t. each column variable
setMethod(".column_grad", "Norm1", function(object, value) {
  rows <- size(object@args[[1]])
  D_null <- Matrix(0, nrow = rows, ncol = 1, sparse = TRUE)
  D_null <- D_null + (value > 0)
  D_null <- D_null - (value < 0)
  D_null   # TODO: Check this is same as ravel and transpose command in CVXPY.
})

#'
#' The NormInf class.
#'
#' This class represents the infinity-norm.
#'
#' @name NormInf-class
#' @aliases NormInf
#' @rdname NormInf-class
.NormInf <- setClass("NormInf", contains = "AxisAtom")

NormInf <- function(x, axis = NA_real_, keepdims = FALSE) { .NormInf(expr = x, axis = axis, keepdims = keepdims) }

#' @param x,object A \linkS4class{NormInf} object.
#' @describeIn NormInf The name and arguments of the atom.
setMethod("name", "NormInf", function(x) {
  paste(class(x), "(", name(x@args[[1]]), ")", sep = "")
})

#' @describeIn NormInf Returns the infinity norm of \code{x}.
setMethod("to_numeric", "NormInf", function(object, values) {
  if(is.na(object@axis))
    # base::norm(values[[1]], type = "I")
    max(abs(values[[1]]))
  else
    # apply_with_keepdims(values[[1]], function(x) { norm(as.matrix(x), type = "I") }, axis = object@axis, keepdims = object@keepdims)
    apply_with_keepdims(values[[1]], function(x) { max(abs(x)) }, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn NormInf Does the atom handle complex numbers?
setMethod("allow_complex", "NormInf", function(object) { TRUE })

#' @describeIn NormInf The atom is always positive.
setMethod("sign_from_args", "NormInf", function(object) { c(TRUE, FALSE) })

#' @describeIn NormInf The atom is convex.
setMethod("is_atom_convex", "NormInf", function(object) { TRUE })

#' @describeIn NormInf The atom is not concave.
setMethod("is_atom_concave", "NormInf", function(object) { FALSE })

#' @describeIn NormInf Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "NormInf", function(object) { TRUE })

#' @describeIn NormInf Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "NormInf", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn NormInf Is the composition weakly increasing in argument \code{idx}?
setMethod("is_incr", "NormInf", function(object, idx) { is_nonneg(object@args[[1]]) })

#' @param idx An index into the atom.
#' @describeIn NormInf Is the composition weakly decreasing in argument \code{idx}?
setMethod("is_decr", "NormInf", function(object, idx) { is_nonpos(object@args[[1]]) })

#' @describeIn NormInf Is the atom piecewise linear?
setMethod("is_pwl", "NormInf", function(object) { is_pwl(object@args[[1]]) })

#' @describeIn NormInf Returns the axis.
setMethod("get_data", "NormInf", function(object) { list(object@axis) })

#' @describeIn NormInf Returns constraints describing the domain of the node
setMethod(".domain", "NormInf", function(object) { list() })

#' @param values A list of numeric values for the arguments
#' @describeIn NormInf Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "NormInf", function(object, values) { .axis_grad(object, values) })

#' @param value A numeric value
#' @describeIn NormInf Gives the (sub/super)gradient of the atom w.r.t. each column variable
setMethod(".column_grad", "NormInf", function(object, value) { stop("Unimplemented") })   # TODO: Implement this! })

#'
#' The NormNuc class.
#'
#' The nuclear norm, i.e. sum of the singular values of a matrix.
#'
#' @slot A An \linkS4class{Expression} or numeric matrix.
#' @name NormNuc-class
#' @aliases NormNuc
#' @rdname NormNuc-class
.NormNuc <- setClass("NormNuc", representation(A = "ConstValORExpr"), contains = "Atom")

#' @param A An \linkS4class{Expression} or numeric matrix.
#' @rdname NormNuc-class
NormNuc <- function(A) { .NormNuc(A = A) }

setMethod("initialize", "NormNuc", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., atom_args = list(.Object@A))
})

#' @param object A \linkS4class{NormNuc} object.
#' @param values A list of arguments to the atom.
#' @describeIn NormNuc The nuclear norm (i.e., the sum of the singular values) of \code{A}.
setMethod("to_numeric", "NormNuc", function(object, values) {
  # Returns the nuclear norm (i.e. the sum of the singular values) of A
  sum(svd(values[[1]])$d)
})

#' @describeIn NormNuc Does the atom handle complex numbers?
setMethod("allow_complex", "NormNuc", function(object) { TRUE })

#' @describeIn NormNuc The atom is a scalar.
setMethod("dim_from_args", "NormNuc", function(object) { c(1,1) })

#' @describeIn NormNuc The atom is positive.
setMethod("sign_from_args",  "NormNuc", function(object) { c(TRUE, FALSE) })

#' @describeIn NormNuc The atom is convex.
setMethod("is_atom_convex", "NormNuc", function(object) { TRUE })

#' @describeIn NormNuc The atom is not concave.
setMethod("is_atom_concave", "NormNuc", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn NormNuc The atom is not monotonic in any argument.
setMethod("is_incr", "NormNuc", function(object, idx) { FALSE })

#' @param idx An index into the atom.
#' @describeIn NormNuc The atom is not monotonic in any argument.
setMethod("is_decr", "NormNuc", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn NormNuc Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "NormNuc", function(object, values) {
  # Grad: UV^T
  s <- svd(values[[1]])
  D <- s$u %*% t(s$v)
  list(Matrix(as.vector(D), sparse = TRUE))
})

#'
#' The OneMinusPos class.
#'
#' This class represents the difference \eqn{1 - x} with domain \eqn{\{x : 0 < x < 1}\}
#'
#' @slot x An \linkS4class{Expression} or numeric matrix.
#' @name OneMinusPos-class
#' @aliases OneMinusPos
#' @rdname OneMinusPos-class
.OneMinusPos <- setClass("OneMinusPos", representation(x = "ConstValORExpr", .ones = "ConstVal"), prototype(.ones = NA_real_), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric matrix.
#' @rdname OneMinusPos-class
OneMinusPos <- function(x) { .OneMinusPos(x = x) }

setMethod("initialize", "OneMinusPos", function(.Object, ..., x) {
  .Object@x <- x
  .Object@.ones <- matrix(1, nrow = nrow(x), ncol = ncol(x))
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@x))
  .Object@args[[1]] <- x
  .Object
})

#' @describeIn OneMinusPos The name and arguments of the atom.
setMethod("name", "OneMinusPos", function(x) { paste(class(x), x@args[[1]]) })

#' @param object A \linkS4class{OneMinusPos} object.
#' @param values A list of arguments to the atom.
#' @describeIn OneMinusPos Returns one minus the value.
setMethod("to_numeric", "OneMinusPos", function(object, values) { object@.ones - values[[1]] })

#' @describeIn OneMinusPos The dimensions of the atom.
setMethod("dim_from_args", "OneMinusPos", function(object) { dim(object@args[[1]]) })

#' @describeIn OneMinusPos Returns the sign (is positive, is negative) of the atom.
setMethod("sign_from_args", "OneMinusPos", function(object) { c(TRUE, FALSE) })

#' @describeIn OneMinusPos Is the atom convex?
setMethod("is_atom_convex", "OneMinusPos", function(object) { FALSE })

#' @describeIn OneMinusPos Is the atom concave?
setMethod("is_atom_concave", "OneMinusPos", function(object) { FALSE })

#' @describeIn OneMinusPos Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "OneMinusPos", function(object) { FALSE })

#' @describeIn OneMinusPos Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "OneMinusPos", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn OneMinusPos Is the atom weakly increasing in the argument \code{idx}?
setMethod("is_incr", "OneMinusPos", function(object, idx) { FALSE })

#' @param idx An index into the atom.
#' @describeIn OneMinusPos Is the atom weakly decreasing in the argument \code{idx}?
setMethod("is_decr", "OneMinusPos", function(object, idx) { TRUE })

#' @param values A list of numeric values for the arguments
#' @describeIn OneMinusPos Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "OneMinusPos", function(object, values) { Matrix(-object@.ones, sparse = TRUE) })

#'
#' The DiffPos atom.
#'
#' The difference between expressions, \eqn{x - y}, where \eqn{x > y > 0}.
#'
#' @param x An \linkS4class{Expression}
#' @param y An \linkS4class{Expression}
#' @return The difference x - y with domain {x,y: x > y > 0}.
DiffPos <- function(x, y) {
  x * OneMinusPos(y/x)
}

#'
#' The PfEigenvalue class.
#'
#' This class represents the Perron-Frobenius eigenvalue of a positive matrix.
#'
#' @slot X An \linkS4class{Expression} or numeric matrix.
#' @name PfEigenvalue-class
#' @aliases PfEigenvalue
#' @rdname PfEigenvalue-class
.PfEigenvalue <- setClass("PfEigenvalue", representation(X = "ConstValORExpr"),
                          validity = function(object) {
                            if(length(dim(object@X)) != 2 || nrow(object@X) != ncol(object@X))
                              stop("[PfEigenvalue: X] The argument X must be a square matrix")
                            return(TRUE)
                          }, contains = "Atom")

#' @param X An \linkS4class{Expression} or numeric matrix.
#' @rdname PfEigenvalue-class
PfEigenvalue <- function(X) { .PfEigenvalue(X = X) }

setMethod("initialize", "PfEigenvalue", function(.Object, ..., X = X) {
  .Object@X <- X
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@X))
  .Object@args[[1]] <- X
  .Object
})

#' @param x,object A \linkS4class{PfEigenvalue} object.
#' @describeIn PfEigenvalue The name and arguments of the atom.
setMethod("name", "PfEigenvalue", function(x) { paste(class(x), x@args[[1]]) })

#' @param values A list of arguments to the atom.
#' @describeIn PfEigenvalue Returns the Perron-Frobenius eigenvalue of \code{X}.
setMethod("to_numeric", "PfEigenvalue", function(object, values) {
  eig <- eigen(values[[1]], only.values = TRUE)
  max(abs(eig$values))
})

#' @describeIn PfEigenvalue The dimensions of the atom.
setMethod("dim_from_args", "PfEigenvalue", function(object) { c(1,1) })

#' @describeIn PfEigenvalue Returns the sign (is positive, is negative) of the atom.
setMethod("sign_from_args", "PfEigenvalue", function(object) { c(TRUE, FALSE) })

#' @describeIn PfEigenvalue Is the atom convex?
setMethod("is_atom_convex", "PfEigenvalue", function(object) { FALSE })

#' @describeIn PfEigenvalue Is the atom concave?
setMethod("is_atom_concave", "PfEigenvalue", function(object) { FALSE })

#' @describeIn PfEigenvalue Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "PfEigenvalue", function(object) { TRUE })

#' @describeIn PfEigenvalue Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "PfEigenvalue", function(object) { FALSE })

# TODO: Figure out monotonicity.
#' @param idx An index into the atom.
#' @describeIn PfEigenvalue Is the atom weakly increasing in the argument \code{idx}?
setMethod("is_incr", "PfEigenvalue", function(object, idx) { FALSE })

#' @param idx An index into the atom.
#' @describeIn PfEigenvalue Is the atom weakly decreasing in the argument \code{idx}?
setMethod("is_decr", "PfEigenvalue", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn PfEigenvalue Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "PfEigenvalue", function(object, values) { NA_real_ })

#'
#' The ProdEntries class.
#'
#' The product of the entries in an expression.
#'
#' @slot expr An \linkS4class{Expression} representing a vector or matrix.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @name ProdEntries-class
#' @aliases ProdEntries
#' @rdname ProdEntries-class
.ProdEntries <- setClass("ProdEntries", contains = "AxisAtom")

#' @param ... \linkS4class{Expression} objects, vectors, or matrices.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @rdname ProdEntries-class
ProdEntries <- function(..., axis = NA_real_, keepdims = FALSE) {
  exprs <- list(...)
  if(length(exprs) == 0)
    stop("Must provide at least one expression")
  else if(length(exprs) == 1)
    .ProdEntries(expr = exprs[[1]], axis = axis, keepdims = keepdims)
  else
    .ProdEntries(expr = do.call("HStack", exprs))
}

#' @param object A \linkS4class{ProdEntries} object.
#' @param values A list of values to take the product of.
#' @describeIn ProdEntries The product of all the entries.
setMethod("to_numeric", "ProdEntries", function(object, values) {
  apply_with_keepdims(values[[1]], prod, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn ProdEntries Returns the sign (is positive, is negative) of the atom.
setMethod("sign_from_args", "ProdEntries", function(object) {
  if(is_nonneg(object@args[[1]]))
    c(TRUE, FALSE)
  else
    c(FALSE, FALSE)
})

#' @describeIn ProdEntries Is the atom convex?
setMethod("is_atom_convex", "ProdEntries", function(object) { FALSE })

#' @describeIn ProdEntries Is the atom concave?
setMethod("is_atom_concave", "ProdEntries", function(object) { FALSE })

#' @describeIn ProdEntries Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "ProdEntries", function(object) { TRUE })

#' @describeIn ProdEntries is the atom log-log concave?
setMethod("is_atom_log_log_concave", "ProdEntries", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn ProdEntries Is the atom weakly increasing in the argument \code{idx}?
setMethod("is_incr", "ProdEntries", function(object, idx) { is_nonneg(object@args[[1]]) })

#' @param idx An index into the atom.
#' @describeIn ProdEntries Is the atom weakly decreasing in the argument \code{idx}?
setMethod("is_decr", "ProdEntries", function(object, idx) { FALSE })

#' @param value A numeric value.
#' @describeIn ProdEntries Gives the (sub/super)gradient of the atom w.r.t. each column variable
setMethod(".column_grad", "ProdEntries", function(object, value) { prod(value)/value })

#' @param values A list of numeric values for the arguments
#' @describeIn ProdEntries Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "ProdEntries", function(object, values) { .axis_grad(object, values) })

#'
#' The QuadForm class.
#'
#' This class represents the quadratic form \eqn{x^T P x}
#'
#' @slot x An \linkS4class{Expression} or numeric vector.
#' @slot P An \linkS4class{Expression}, numeric matrix, or vector.
#' @name QuadForm-class
#' @aliases QuadForm
#' @rdname QuadForm-class
.QuadForm <- setClass("QuadForm", representation(x = "ConstValORExpr", P = "ConstValORExpr"), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric vector.
#' @param P An \linkS4class{Expression}, numeric matrix, or vector.
#' @rdname QuadForm-class
QuadForm <- function(x, P) { .QuadForm(x = x, P = P) }

setMethod("initialize", "QuadForm", function(.Object, ..., x, P) {
  .Object@x <- x
  .Object@P <- P
  callNextMethod(.Object, ..., atom_args = list(.Object@x, .Object@P))
})

#' @describeIn QuadForm The name and arguments of the atom.
setMethod("name", "QuadForm", function(x) {
  paste(class(x), "(", x@args[[1]], ", ", x@args[[2]], ")", sep = "")
})

#' @param object A \linkS4class{QuadForm} object.
#' @describeIn QuadForm Does the atom handle complex numbers?
setMethod("allow_complex", "QuadForm", function(object) { TRUE })

#' @param values A list of numeric values for the arguments
#' @describeIn QuadForm Returns the quadratic form.
setMethod("to_numeric", "QuadForm", function(object, values) {
  prod <- values[[2]] %*% values[[1]]
  if(is_complex(object@args[[1]]))
    return(t(Conj(values[[1]])) %*% prod)
  else
    return(t(values[[1]]) %*% prod)
})

#' @describeIn QuadForm Checks the dimensions of the arguments.
setMethod("validate_args", "QuadForm", function(object) {
  callNextMethod()
  n <- nrow(object@args[[2]])
  x_dim <- dim(object@args[[1]])
  # if(ncol(object@args[[2]]) != n || !(dim(object@args[[1]]) %in% list(c(n, 1), c(n, NA_real_))))
  if(ncol(object@args[[2]]) != n || !(length(x_dim) == 2 && all(x_dim == c(n,1))))
    stop("Invalid dimensions for arguments.")
})

#' @describeIn QuadForm Returns the sign (is positive, is negative) of the atom.
setMethod("sign_from_args", "QuadForm", function(object) { c(is_atom_convex(object), is_atom_concave(object)) })

#' @describeIn QuadForm The dimensions of the atom.
setMethod("dim_from_args", "QuadForm", function(object) {
  # if(ndim(object@args[[1]]) == 0)
  #  c()
  # else
  #  c(1,1)
  c(1,1)
})

#' @describeIn QuadForm Is the atom convex?
setMethod("is_atom_convex", "QuadForm", function(object) { is_psd(object@args[[2]]) })

#' @describeIn QuadForm Is the atom concave?
setMethod("is_atom_concave", "QuadForm", function(object) { is_nsd(object@args[[2]]) })

#' @describeIn QuadForm Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "QuadForm", function(object) { TRUE })

#' @describeIn QuadForm Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "QuadForm", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn QuadForm Is the atom weakly increasing in the argument \code{idx}?
setMethod("is_incr", "QuadForm", function(object, idx) {
  (is_nonneg(object@args[[1]]) && is_nonneg(object@args[[2]])) ||
    (is_nonpos(object@args[[1]]) && is_nonneg(object@args[[2]]))
})

#' @param idx An index into the atom.
#' @describeIn QuadForm Is the atom weakly decreasing in the argument \code{idx}?
setMethod("is_decr", "QuadForm", function(object, idx) {
  (is_nonneg(object@args[[1]]) && is_nonpos(object@args[[2]])) ||
    (is_nonpos(object@args[[1]]) && is_nonpos(object@args[[2]]))
})

#' @describeIn QuadForm Is the atom quadratic?
setMethod("is_quadratic", "QuadForm", function(object) { TRUE })

#' @describeIn QuadForm Is the atom piecewise linear?
setMethod("is_pwl", "QuadForm", function(object) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn QuadForm Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "QuadForm", function(object, values) {
  x <- values[[1]]
  P <- values[[2]]
  D <- 2*P %*% t(x)
  Matrix(as.vector(t(D)), sparse = TRUE)
})

#'
#' The SymbolicQuadForm class.
#'
#' @slot x An \linkS4class{Expression} or numeric vector.
#' @slot P An \linkS4class{Expression}, numeric matrix, or vector.
#' @slot original_expression The original \linkS4class{Expression}.
#' @name SymbolicQuadForm-class
#' @aliases SymbolicQuadForm
#' @rdname SymbolicQuadForm-class
.SymbolicQuadForm <- setClass("SymbolicQuadForm", representation(x = "ConstValORExpr", P = "ConstValORExpr", original_expression = "Expression"), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric vector.
#' @param P An \linkS4class{Expression}, numeric matrix, or vector.
#' @param expr The original \linkS4class{Expression}.
#' @rdname SymbolicQuadForm-class
SymbolicQuadForm <- function(x, P, expr) { .SymbolicQuadForm(x = x, P = P, original_expression = expr) }

setMethod("initialize", "SymbolicQuadForm", function(.Object, ..., x, P, original_expression) {
  .Object@x <- x
  .Object@original_expression <- original_expression
  .Object <- callNextMethod(.Object, ..., atom_args = list(x, P), validate = FALSE)
  .Object@P <- .Object@args[[2]]
  validObject(.Object)
  .Object
})

#' @param object A \linkS4class{SymbolicQuadForm} object.
#' @describeIn SymbolicQuadForm The dimensions of the atom.
setMethod("dim_from_args", "SymbolicQuadForm", function(object) { dim_from_args(object@original_expression) })

#' @describeIn SymbolicQuadForm The sign (is positive, is negative) of the atom.
setMethod("sign_from_args", "SymbolicQuadForm", function(object) { sign_from_args(object@original_expression) })

#' @describeIn SymbolicQuadForm The original expression.
setMethod("get_data", "SymbolicQuadForm", function(object) { list(object@original_expression) })

#' @describeIn SymbolicQuadForm Is the original expression convex?
setMethod("is_atom_convex", "SymbolicQuadForm", function(object) { is_atom_convex(object@original_expression) })

#' @describeIn SymbolicQuadForm Is the original expression concave?
setMethod("is_atom_concave", "SymbolicQuadForm", function(object) { is_atom_concave(object@original_expression) })

#' @param idx An index into the atom.
#' @describeIn SymbolicQuadForm Is the original expression weakly increasing in argument \code{idx}?
setMethod("is_incr", "SymbolicQuadForm", function(object, idx) { is_incr(object@original_expression, idx) })

#' @describeIn SymbolicQuadForm Is the original expression weakly decreasing in argument \code{idx}?
setMethod("is_decr", "SymbolicQuadForm", function(object, idx) { is_decr(object@original_expression, idx) })

#' @describeIn SymbolicQuadForm The atom is quadratic.
setMethod("is_quadratic", "SymbolicQuadForm", function(object) { TRUE })

#' @param values A list of numeric values for the arguments
#' @describeIn SymbolicQuadForm Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "SymbolicQuadForm", function(object, values) { stop("Unimplemented") })

#'
#' Compute a Matrix Decomposition.
#'
#' Compute sgn, scale, M such that \eqn{P = sgn * scale * dot(M, t(M))}.
#'
#' @param P A real symmetric positive or negative (semi)definite input matrix
#' @param cond Cutoff for small eigenvalues. Singular values smaller than rcond * largest_eigenvalue are considered negligible.
#' @param rcond Cutoff for small eigenvalues. Singular values smaller than rcond * largest_eigenvalue are considered negligible.
#' @return A list consisting of induced matrix 2-norm of P and a rectangular matrix such that P = scale * (dot(M1, t(M1)) - dot(M2, t(M2)))
.decomp_quad <- function(P, cond = NA, rcond = NA) {
  eig <- eigen(P, only.values = FALSE)
  w <- eig$values
  V <- eig$vectors

  if(!is.na(rcond))
    cond <- rcond
  if(cond == -1 || is.na(cond))
    cond <- 1e6 * .Machine$double.eps   # All real numbers are stored as double precision in R

  scale <- max(abs(w))
  if(scale < cond)
    return(list(scale = 0, M1 = V[,FALSE], M2 = V[,FALSE]))
  w_scaled <- w / scale
  maskp <- w_scaled > cond
  maskn <- w_scaled < -cond

  # TODO: Allow indefinite QuadForm
  if(any(maskp) && any(maskn))
    warning("Forming a non-convex expression QuadForm(x, indefinite)")

  if(sum(maskp) <= 1)
    M1 <- as.matrix(V[,maskp] * sqrt(w_scaled[maskp]))
  else
    M1 <- V[,maskp] %*% diag(sqrt(w_scaled[maskp]))

  if(sum(maskn) <= 1)
    M2 <- as.matrix(V[,maskn]) * sqrt(-w_scaled[maskn])
  else
    M2 <- V[,maskn] %*% diag(sqrt(-w_scaled[maskn]))
  list(scale = scale, M1 = M1, M2 = M2)
}

#'
#' The QuadOverLin class.
#'
#' This class represents the sum of squared entries in X divided by a scalar y, \eqn{\sum_{i,j} X_{i,j}^2/y}.
#'
#' @slot x An \linkS4class{Expression} or numeric matrix.
#' @slot y A scalar \linkS4class{Expression} or numeric constant.
#' @name QuadOverLin-class
#' @aliases QuadOverLin
#' @rdname QuadOverLin-class
.QuadOverLin <- setClass("QuadOverLin", representation(x = "ConstValORExpr", y = "ConstValORExpr"), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric matrix.
#' @param y A scalar \linkS4class{Expression} or numeric constant.
#' @rdname QuadOverLin-class
QuadOverLin <- function(x, y) { .QuadOverLin(x = x, y = y) }

setMethod("initialize", "QuadOverLin", function(.Object, ..., x, y) {
  .Object@x <- x
  .Object@y <- y
  callNextMethod(.Object, ..., atom_args = list(.Object@x, .Object@y))
})

#' @describeIn QuadOverLin Does the atom handle complex numbers?
setMethod("allow_complex", "QuadOverLin", function(object) { TRUE })

#' @param object A \linkS4class{QuadOverLin} object.
#' @param values A list of arguments to the atom.
#' @describeIn QuadOverLin The sum of the entries of \code{x} squared over \code{y}.
setMethod("to_numeric", "QuadOverLin", function(object, values) { sum(Mod(values[[1]])^2) / values[[2]] })

#' @describeIn QuadOverLin Check the dimensions of the arguments.
setMethod("validate_args",   "QuadOverLin", function(object) {
  if(!is_scalar(object@args[[2]]))
    stop("The second argument to QuadOverLin must be a scalar.")
  if(is_complex(object@args[[2]]))
    stop("The second argument to QuadOverLin cannot be complex.")
  callNextMethod()
})

#' @describeIn QuadOverLin The atom is a scalar.
setMethod("dim_from_args", "QuadOverLin", function(object) { c(1,1) })

#' @describeIn QuadOverLin The atom is positive.
setMethod("sign_from_args",  "QuadOverLin", function(object) { c(TRUE, FALSE) })

#' @describeIn QuadOverLin The atom is convex.
setMethod("is_atom_convex", "QuadOverLin", function(object) { TRUE })

#' @describeIn QuadOverLin The atom is not concave.
setMethod("is_atom_concave", "QuadOverLin", function(object) { FALSE })

#' @describeIn QuadOverLin Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "QuadOverLin", function(object) { TRUE })

#' @describeIn QuadOverLin Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "QuadOverLin", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn QuadOverLin A logical value indicating whether the atom is weakly increasing in argument \code{idx}.
setMethod("is_incr", "QuadOverLin", function(object, idx) { (idx == 1) && is_nonneg(object@args[[idx]]) })

#' @describeIn QuadOverLin A logical value indicating whether the atom is weakly decreasing in argument \code{idx}.
setMethod("is_decr", "QuadOverLin", function(object, idx) { ((idx == 1) && is_nonpos(object@args[[idx]])) || (idx == 2) })

#' @describeIn QuadOverLin Quadratic if \code{x} is affine and \code{y} is constant.
setMethod("is_quadratic", "QuadOverLin", function(object) { is_affine(object@args[[1]]) && is_constant(object@args[[2]]) })

#' @describeIn QuadOverLin Quadratic of piecewise affine if \code{x} is piecewise linear and \code{y} is constant.
setMethod("is_qpwa", "QuadOverLin", function(object) { is_pwl(object@args[[1]]) && is_constant(object@args[[2]]) })

#' @describeIn QuadOverLin Returns constraints describing the domain of the node
setMethod(".domain", "QuadOverLin", function(object) { list(object@args[[2]] >= 0) })

#' @param values A list of numeric values for the arguments
#' @describeIn QuadOverLin Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "QuadOverLin", function(object, values) {
  X <- values[[1]]
  y <- as.vector(values[[2]])
  if(y <= 0)
    return(list(NA_real_, NA_real_))
  else {
    # DX = 2X/y, Dy = -||X||^2_2/y^2
    Dy <- -sum(Mod(X)^2)/y^2
    Dy <- Matrix(Dy, sparse = TRUE)
    DX <- 2.0*X/y
    DX <- Matrix(as.vector(t(DX)), sparse = TRUE)
    return(list(DX, Dy))
  }
})

#'
#' The SigmaMax class.
#'
#' The maximum singular value of a matrix.
#'
#' @slot A An \linkS4class{Expression} or numeric matrix.
#' @name SigmaMax-class
#' @aliases SigmaMax
#' @rdname SigmaMax-class
.SigmaMax <- setClass("SigmaMax", representation(A = "ConstValORExpr"), contains = "Atom")

#' @param A An \linkS4class{Expression} or matrix.
#' @rdname SigmaMax-class
SigmaMax <- function(A = A) { .SigmaMax(A = A) }

setMethod("initialize", "SigmaMax", function(.Object, ..., A) {
  .Object@A <- A
  callNextMethod(.Object, ..., atom_args = list(.Object@A))
})

#' @param object A \linkS4class{SigmaMax} object.
#' @param values A list of arguments to the atom.
#' @describeIn SigmaMax The largest singular value of \code{A}.
setMethod("to_numeric", "SigmaMax", function(object, values) { base::norm(values[[1]], type = "2") })

#' @describeIn SigmaMax Does the atom handle complex numbers?
setMethod("allow_complex", "SigmaMax", function(object) { TRUE })

#' @describeIn SigmaMax The atom is a scalar.
setMethod("dim_from_args", "SigmaMax", function(object) { c(1,1) })

#' @describeIn SigmaMax The atom is positive.
setMethod("sign_from_args",  "SigmaMax", function(object) { c(TRUE, FALSE) })

#' @describeIn SigmaMax The atom is convex.
setMethod("is_atom_convex", "SigmaMax", function(object) { TRUE })

#' @describeIn SigmaMax The atom is concave.
setMethod("is_atom_concave", "SigmaMax", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn SigmaMax The atom is not monotonic in any argument.
setMethod("is_incr", "SigmaMax", function(object, idx) { FALSE })

#' @describeIn SigmaMax The atom is not monotonic in any argument.
setMethod("is_decr", "SigmaMax", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn SigmaMax Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "SigmaMax", function(object, values) {
  # Grad: U diag(e_1) t(V)
  s <- svd(values[[1]])
  ds <- rep(0, length(s$d))
  ds[1] <- 1
  D <- s$u %*% diag(ds) %*% t(s$v)
  list(Matrix(as.vector(D), sparse = TRUE))
})

#'
#' The SumLargest class.
#'
#' The sum of the largest k values of a matrix.
#'
#' @slot x An \linkS4class{Expression} or numeric matrix.
#' @slot k The number of largest values to sum over.
#' @name SumLargest-class
#' @aliases SumLargest
#' @rdname SumLargest-class
.SumLargest <- setClass("SumLargest", representation(x = "ConstValORExpr", k = "numeric"), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric matrix.
#' @param k The number of largest values to sum over.
#' @rdname SumLargest-class
SumLargest <- function(x, k) { .SumLargest(x = x, k = k) }

setMethod("initialize", "SumLargest", function(.Object, ..., x, k) {
  .Object@x <- x
  .Object@k <- k
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object A \linkS4class{SumLargest} object.
#' @param values A list of arguments to the atom.
#' @describeIn SumLargest The sum of the \code{k} largest entries of the vector or matrix.
setMethod("to_numeric", "SumLargest", function(object, values) {
  # Return the sum of the k largest entries of the matrix
  value <- as.vector(values[[1]])
  k <- min(object@k, length(value))
  val_sort <- sort(value, decreasing = TRUE)
  sum(val_sort[1:k])
})

#' @describeIn SumLargest Check that \code{k} is a positive integer.
setMethod("validate_args",   "SumLargest", function(object) {
  if(as.integer(object@k) != object@k || object@k <= 0)
    stop("[SumLargest: validation] k must be a positive integer")
  callNextMethod()
})

#' @describeIn SumLargest The atom is a scalar.
setMethod("dim_from_args", "SumLargest", function(object) { c(1,1) })

#' @describeIn SumLargest The sign of the atom.
setMethod("sign_from_args", "SumLargest", function(object) { c(is_nonneg(object@args[[1]]), is_nonpos(object@args[[1]])) })

#' @describeIn SumLargest The atom is convex.
setMethod("is_atom_convex", "SumLargest", function(object) { TRUE })

#' @describeIn SumLargest The atom is not concave.
setMethod("is_atom_concave", "SumLargest", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn SumLargest The atom is weakly increasing in every argument.
setMethod("is_incr", "SumLargest", function(object, idx) { TRUE })

#' @describeIn SumLargest The atom is not weakly decreasing in any argument.
setMethod("is_decr", "SumLargest", function(object, idx) { FALSE })

#' @describeIn SumLargest A list containing \code{k}.
setMethod("get_data", "SumLargest", function(object) { list(object@k) })

#' @param values A list of numeric values for the arguments
#' @describeIn SumLargest Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "SumLargest", function(object, values) {
  # Grad: 1 for each of the k largest indices
  value <- as.vector(t(values[[1]]))
  k <- min(object@k, length(value))
  indices <- order(value, decreasing = TRUE)
  arg_dim <- dim(object@args[[1]])
  D <- matrix(0, nrow = arg_dim[1]*arg_dim[2], ncol = 1)
  D[indices[1:k]] <- 1
  list(Matrix(D, sparse = TRUE))
})

#'
#' The SumSmallest atom.
#'
#' The sum of the smallest k values of a matrix.
#'
#' @param x An \linkS4class{Expression} or numeric matrix.
#' @param k The number of smallest values to sum over.
#' @return Sum of the smlalest k values
SumSmallest <- function(x, k) {
  x <- as.Constant(x)
  -SumLargest(x = -x, k = k)
}

#'
#' The SumSquares atom.
#'
#' The sum of the squares of the entries.
#'
#' @param expr An \linkS4class{Expression} or numeric matrix.
#' @return Sum of the squares of the entries in the expression.
SumSquares <- function(expr) { QuadOverLin(x = expr, y = 1) }

#'
#' The TotalVariation atom.
#'
#' The total variation of a vector, matrix, or list of matrices.
#' Uses L1 norm of discrete gradients for vectors and L2 norm of discrete gradients for matrices.
#'
#' @param value An \linkS4class{Expression} representing the value to take the total variation of.
#' @param ... Additional matrices extending the third dimension of value.
#' @return An expression representing the total variation.
TotalVariation <- function(value, ...) {
  value <- as.Constant(value)
  if(ndim(value) == 0 || (nrow(value) == 1 && ncol(value) == 1))
    stop("TotalVariation cannot take a scalar argument")
  else if(ndim(value) == 1 || nrow(value) == 1 || ncol(value) == 1) {  # L1 norm for vectors
    if(nrow(value) == 1)
      value <- t(value)
    Norm(value[2:nrow(value),1] - value[1:(nrow(value)-1),1], 1)
  } else {   # L2 norm for matrices
    val_dim <- dim(value)
    rows <- val_dim[1]
    cols <- val_dim[2]
    args <- lapply(list(...), as.Constant)
    values <- c(list(value), args)

    diffs <- list()
    for(mat in values) {
      if(rows > 1 && cols > 1) {
        diffs <- c(diffs, list(mat[1:(rows-1), 2:cols] - mat[1:(rows-1), 1:(cols-1)],
                               mat[2:rows, 1:(cols-1)] - mat[1:(rows-1), 1:(cols-1)]))
      } else {
        diffs <- c(diffs, list(matrix(0, nrow = rows-1, ncol = cols-1),
                               matrix(0, nrow = rows-1, ncol = cols-1)))
      }
    }
    len <- nrow(diffs[[1]]) * ncol(diffs[[2]])
    stacked <- .VStack(atom_args = lapply(diffs, function(diff) { Reshape(diff, c(1, len)) }))
    SumEntries(Norm(stacked, p = 2, axis = 2))
  }
}
