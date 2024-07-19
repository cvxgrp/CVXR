#'
#' The Leaf class.
#'
#' This class represents a leaf node, i.e. a Variable, Constant, or Parameter.
#'
#' @slot id (Internal) A unique integer identification number used internally.
#' @slot dim The dimensions of the leaf.
#' @slot value The numeric value of the leaf.
#' @slot nonneg Is the leaf constrained to be nonnegative?
#' @slot nonpos Is the leaf constrained to be nonpositive?
#' @slot complex Is the leaf complex valued?
#' @slot imag Is the leaf imaginary?
#' @slot symmetric Is the leaf symmetric?
#' @slot diag Is the leaf diagonal?
#' @slot PSD Is the leaf constrained to be positive semidefinite?
#' @slot NSD Is the leaf constrained to be negative semidefinite?
#' @slot hermitian Is the leaf Hermitian?
#' @slot boolean Is the leaf boolean? Can be \code{TRUE} = entire leaf is boolean, \code{FALSE} = entire leaf is not boolean, or a vector of
#' indices which should be constrained as boolean, where each index is a vector of length exactly equal to the length of \code{dim}.
#' @slot integer Is the leaf integral? The semantics are the same as the \code{boolean} argument.
#' @slot sparsity A matrix representing the fixed sparsity pattern of the leaf.
#' @slot pos Is the leaf positive?
#' @slot neg Is the leaf negative?
#' @name Leaf-class
#' @aliases Leaf
#' @rdname Leaf-class
Leaf <- setClass("Leaf", representation(dim = "NumORNULL", value = "ConstVal", nonneg = "logical", nonpos = "logical",
                                        complex = "logical", imag = "logical", symmetric = "logical", diag = "logical",
                                        PSD = "logical", NSD = "logical", hermitian = "logical", boolean = "NumORLogical", integer = "NumORLogical",
                                        sparsity = "matrix", pos = "logical", neg = "logical",
                                        attributes = "list", boolean_idx = "matrix", integer_idx = "matrix"),
                         prototype(value = NA_real_, nonneg = FALSE, nonpos = FALSE,
                                   complex = FALSE, imag = FALSE, symmetric = FALSE, diag = FALSE,
                                   PSD = FALSE, NSD = FALSE, hermitian = FALSE, boolean = FALSE, integer = FALSE,
                                   sparsity = matrix(0, nrow = 0, ncol = 1), pos = FALSE, neg = FALSE,
                                   attributes = list(), boolean_idx = matrix(0, nrow = 0, ncol = 1), integer_idx = matrix(0, nrow = 0, ncol = 1)), contains = "Expression")

setMethod("initialize", "Leaf", function(.Object, ..., dim, value = NA_real_, nonneg = FALSE, nonpos = FALSE, complex = FALSE, imag = FALSE, symmetric = FALSE, diag = FALSE, PSD = FALSE, NSD = FALSE, hermitian = FALSE, boolean = FALSE, integer = FALSE, sparsity = matrix(0, nrow = 0, ncol = 1), pos = FALSE, neg = FALSE, attributes = list(), boolean_idx = matrix(0, nrow = 0, ncol = 1), integer_idx = matrix(0, nrow = 0, ncol = 1)) {
  if(length(dim) > 2)
    stop("Expressions of dimension greater than 2 are not supported.")

  for(d in dim) {
    if(!intf_is_integer(d) || d <= 0)
      stop("Invalid dimensions ", dim)
  }
  .Object@dim <- as.integer(dim)

  if((PSD || NSD || symmetric || diag || hermitian) && (length(dim) != 2 || dim[1] != dim[2]))
    stop("Invalid dimensions (", paste(dim, collapse = ", "), "). Must be a square matrix.")

  # Construct matrix of boolean/integer-constrained indices.
  if(is.logical(boolean)) {
    if(boolean)
      .Object@boolean_idx <- as.matrix(do.call(expand.grid, lapply(dim, function(k) { seq(k) })))
    else
      .Object@boolean_idx <- matrix(0, nrow = 0, ncol = length(dim))
    bool_attr <- boolean
  } else {
    .Object@boolean_idx <- boolean
    bool_attr <- (nrow(boolean) > 0)
  }

  if(is.logical(integer)) {
    if(integer)
      .Object@integer_idx <- as.matrix(do.call(expand.grid, lapply(dim, function(k) { 1:k })))
    else
      .Object@integer_idx <- matrix(0, nrow = 0, ncol = length(dim))
    int_attr <- integer
  } else {
    .Object@integer_idx <- integer
    int_attr <- (nrow(integer) > 0)
  }

  # Process attributes.
  .Object@attributes <- list(nonneg = nonneg, nonpos = nonpos, pos = pos, neg = neg, complex = complex, imag = imag,
                             symmetric = symmetric, diag = diag, PSD = PSD, NSD = NSD, hermitian = hermitian,
                             boolean = bool_attr, integer = int_attr, sparsity = sparsity)

  # Only one attribute can be TRUE (except boolean and integer).
  attrs <- .Object@attributes
  attrs$sparsity <- prod(dim(attrs$sparsity)) != 0
  true_attr <- sum(unlist(attrs))

  if(bool_attr && int_attr)
    true_attr <- true_attr - 1
  if(true_attr > 1)
    stop("Cannot set more than one special attribute.")

  if(!any(is.na(value)))
    .Object@value <- value
  callNextMethod(.Object, ..., args = list())
})

setMethod("get_attr_str", "Leaf", function(object) {
  # Get a string representing the attributes
  attr_str <- ""
  for(attr in names(object@attributes)) {
    val <- object@attributes[[attr]]
    if(attr != "real" && !is.null(val)) {
      if(nchar(attr_str) == 0)
        attr_str <- sprintf("%s=%s", attr, val)
      else
        attr_str <- paste(attr_str, sprintf("%s=%s", attr, val), sep = ", ")
    }
  }
  attr_str
})

# TODO: Get rid of this and just skip calling copy on Leaf objects.
setMethod("copy", "Leaf", function(object, args = NULL, id_objects = list()) {
  # if("id" %in% names(attributes(object)) && as.character(object@id) %in% names(id_objects))
  if(!is.na(object@id) && as.character(object@id) %in% names(id_objects))
    return(id_objects[[as.character(object@id)]])
  return(object)   # Leaves are not deep copied.
})

#' @param object,x A \linkS4class{Leaf} object.
#' @describeIn Leaf Leaves are not copied.
setMethod("get_data", "Leaf", function(object) { list() })

#' @describeIn Leaf The dimensions of the leaf node.
setMethod("dim", "Leaf", function(x) { x@dim })

#' @describeIn Leaf List of \linkS4class{Variable} objects in the leaf node.
setMethod("variables", "Leaf", function(object) { list() })

#' @describeIn Leaf List of \linkS4class{Parameter} objects in the leaf node.
setMethod("parameters", "Leaf", function(object) { list() })

#' @describeIn Leaf List of \linkS4class{Constant} objects in the leaf node.
setMethod("constants", "Leaf", function(object) { list() })

#' @describeIn Leaf List of \linkS4class{Atom} objects in the leaf node.
setMethod("atoms", "Leaf", function(object) { list() })

#' @describeIn Leaf A logical value indicating whether the leaf node is convex.
setMethod("is_convex", "Leaf", function(object) { TRUE })

#' @describeIn Leaf A logical value indicating whether the leaf node is concave.
setMethod("is_concave", "Leaf", function(object) { TRUE })

#' @describeIn Leaf Is the expression log-log convex?
setMethod("is_log_log_convex", "Leaf", function(object) { is_pos(object) })

#' @describeIn Leaf Is the expression log-log concave?
setMethod("is_log_log_concave", "Leaf", function(object) { is_pos(object) })

#' @describeIn Leaf A logical value indicating whether the leaf node is nonnegative.
setMethod("is_nonneg", "Leaf", function(object) { object@attributes$nonneg || object@attributes$pos || object@attributes$boolean })

#' @describeIn Leaf A logical value indicating whether the leaf node is nonpositive.
setMethod("is_nonpos", "Leaf", function(object) { object@attributes$nonpos || object@attributes$neg })

#' @describeIn Leaf Is the expression positive?
setMethod("is_pos", "Leaf", function(object) { object@attributes$pos })

#' @describeIn Leaf Is the expression negative?
setMethod("is_neg", "Leaf", function(object) { object@attributes$neg })

#' @describeIn Leaf A logical value indicating whether the leaf node is hermitian.
setMethod("is_hermitian", "Leaf", function(object) {
  (is_real(object) && is_symmetric(object)) || object@attributes$hermitian || is_psd(object) || is_nsd(object)
})

#' @describeIn Leaf A logical value indicating whether the leaf node is symmetric.
setMethod("is_symmetric", "Leaf", function(object) {
  is_scalar(object) || any(sapply(c("diag", "symmetric", "PSD", "NSD"), function(key) { object@attributes[[key]] }))
})

#' @describeIn Leaf A logical value indicating whether the leaf node is imaginary.
setMethod("is_imag", "Leaf", function(object) { object@attributes$imag })

#' @describeIn Leaf A logical value indicating whether the leaf node is complex.
setMethod("is_complex", "Leaf", function(object) {
  object@attributes$complex || is_imag(object) || object@attributes$hermitian
})

#' @describeIn Leaf A list of constraints describing the closure of the region where the leaf node is finite. Default is the full domain.
setMethod("domain", "Leaf", function(object) {
  domain <- list()
  if(object@attributes$nonneg || object@attributes$pos)
    domain <- c(domain, object >= 0)
  else if(object@attributes$nonpos || object@attributes$neg)
    domain <- c(domain, object <= 0)
  else if(object@attributes$PSD)
    domain <- c(domain, object %>>% 0)
  else if(object@attributes$NSD)
    domain <- c(domain, object %<<% 0)
  return(domain)
})

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn Leaf Project value onto the attribute set of the leaf.
setMethod("project", "Leaf", function(object, value) {
  # Only one attribute can be active at once (besides real, nonpos/nonneg, and bool/int).
  if(!is_complex(object))
    value <- Re(value)

  if(object@attributes$nonpos && object@attributes$nonneg)
    return(0*value)
  else if(object@attributes$nonpos || object@attributes$neg)
    return(pmin(value, 0))
  else if(object@attributes$nonneg || object@attributes$pos)
    return(pmax(value, 0))
  else if(object@attributes$imag)
    return(Im(value)*1i)
  else if(object@attributes$complex)
    return(as.complex(value))
  else if(object@attributes$boolean)
    # TODO: Respect the boolean indices.
    return(round(pmax(pmin(value, 1), 0)))
  else if(object@attributes$integer)
    # TODO: Respect the integer indices. Also, variable may be integer in some indices and boolean in others.
    return(round(value))
  else if(object@attributes$diag) {
    val <- diag(value)
    return(sparseMatrix(i = 1:length(val), j = 1:length(val), x = val))
  } else if(object@attributes$hermitian)
    return((value + t(Conj(value)))/2)
  else if(any(sapply(c("symmetric", "PSD", "NSD"), function(key) { object@attributes[[key]] }))) {
    value <- value + t(value)
    value <- value/2
    if(object@attributes$symmetric)
      return(value)

    wV <- eigen(value, symmetric = TRUE, only.values = FALSE)
    w <- wV$values
    V <- wV$vectors

    if(object@attributes$PSD) {
      bad <- w < 0
      if(!any(bad))
        return(value)
      w[bad] <- 0
    } else {   # NSD
      bad <- w > 0
      if(!any(bad))
        return(value)
      w[bad] <- 0
    }
    return((V %*% diag(w)) %*% t(V))
  } else
    return(value)
})
#' @describeIn Leaf Get the value of the leaf.
setMethod("value", "Leaf", function(object) { object@value })

#' @describeIn Leaf Set the value of the leaf.
setReplaceMethod("value", "Leaf", function(object, value) {
  object@value <- validate_val(object, value)
  return(object)
})

#' @describeIn Leaf Project and assign a value to the leaf.
setMethod("project_and_assign", "Leaf", function(object, value) {
  object@value <- project(object, value)
  return(object)
})

#' @param val The assigned value.
#' @describeIn Leaf Check that \code{val} satisfies symbolic attributes of leaf.
setMethod("validate_val", "Leaf", function(object, val) {
  if(!any(is.na(val))) {
    # Convert val to numeric vector/matrix or sparse matrix.
    val <- intf_convert(val)
    if(any(intf_dim(val) != dim(object)))
      stop("Invalid dimensions (", paste(intf_dim(val), collapse = ","), ") for value")
    projection <- project(object, val)   # Might be an R vector/matrix or sparse Matrix object.
    delta <- abs(val - projection)

    if(is(delta, "sparseMatrix")) {
      # Based on current implementation of project, it is not possible for this Leaf to be PSD/NSD and a sparse matrix.
      close_enough <- all(abs(delta@x) <= SPARSE_PROJECTION_TOL)   # Only check for near-equality on nonzero values.
    } else {
      delta <- as.matrix(delta)
      # Need to measure residual in a canonical way.
      if(object@attributes$PSD || object@attributes$NSD) {
        # For PSD/NSD Leafs, we use the largest singular value norm.
        close_enough <- norm(delta, type = "2") <= PSD_NSD_PROJECTION_TOL
      } else {
        # For all other Leafs, we use the infinity norm on the vectorized Leaf.
        close_enough <- all(abs(delta) <= GENERAL_PROJECTION_TOL)
      }
    }

    if(!close_enough) {
      if(object@attributes$nonneg)
        attr_str <- "nonnegative"
      else if(object@attributes$pos)
        attr_str <- "positive"
      else if(object@attributes$nonpos)
        attr_str <- "nonpositive"
      else if(object@attributes$neg)
        attr_str <- "negative"
      else if(object@attributes$diag)
        attr_str <- "diagonal"
      else if(object@attributes$PSD)
        attr_str <- "positive semidefinite"
      else if(object@attributes$NSD)
        attr_str <- "negative semidefinite"
      else if(object@attributes$imag)
        attr_str <- "imaginary"
      else {
        attr_str <- names(object@attributes)[unlist(object@attributes) == 1]
        attr_str <- c(attr_str, "real")[1]
      }
      stop("Value must be ", attr_str)
    }
  }
  return(val)
})

#' @describeIn Leaf A logical value indicating whether the leaf node is a positive semidefinite matrix.
setMethod("is_psd", "Leaf", function(object) { object@attributes$PSD })

#' @describeIn Leaf A logical value indicating whether the leaf node is a negative semidefinite matrix.
setMethod("is_nsd", "Leaf", function(object) { object@attributes$NSD })

#' @describeIn Leaf Leaf nodes are always quadratic.
setMethod("is_quadratic", "Leaf", function(object) { TRUE })

#' @describeIn Leaf Leaf nodes are not quadratic terms.
setMethod("has_quadratic_term", "Leaf", function(object) { FALSE })

#' @describeIn Leaf Leaf nodes are always piecewise linear.
setMethod("is_pwl", "Leaf", function(object) { TRUE })

#' @describeIn Leaf Is the expression a disciplined parametrized expression?
setMethod("is_dpp", "Leaf", function(object, context = "dcp") { TRUE })

#' @describeIn Leaf A list of the atoms in the expression.
setMethod("atoms", "Leaf", function(object) { list() })
