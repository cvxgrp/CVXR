## LinOp class shadowing CPP class
CVXcanon.LinOp <- function(type = NULL, size = NULL, args = NULL, data = NULL) {
  operatorType <- c(
    "VARIABLE",
    "PROMOTE",
    "MUL_EXPR",
    "RMUL_EXPR",
    "MUL_ELEM",
    "DIV",
    "SUM",
    "NEG",
    "INDEX",
    "TRANSPOSE",
    "SUM_ENTRIES",
    "TRACE",
    "RESHAPE_EXPR",
    "DIAG_VEC",
    "DIAG_MAT",
    "UPPER_TRI",
    "CONV",
    "HSTACK",
    "VSTACK",
    "SCALAR_CONST",
    "DENSE_CONST",
    "SPARSE_CONST",
    "NO_OP",
    "KRON")
  ## Create a new LinOp on the C side
  ptr <- .Call("_CVXR_LinOp__new", PACKAGE = "CVXR")

  getXPtr <- function() ptr

  ## Associate args on R side with the args on the C side.
  if (!is.null(type)) {
    set_type(type)
  }

  if (!is.null(size)) {
    set_size(size)
  }

  ## make a direct call for efficiency rather than calling self$args_push_back(x)
  for (x in args) .Call("_CVXR_LinOp__args_push_back", ptr, x$getXPtr(), PACKAGE = "CVXR")

  if (!is.null(data)) {
    set_dense_data(data)
  }
  
  set_sparse <- function(value) .Call("_CVXR_LinOp__set_sparse",  ptr, value, PACKAGE = "CVXR")
  get_sparse <- function() .Call("_CVXR_LinOp__get_sparse", ptr, PACKAGE = "CVXR")
  set_sparse_data <- function(value) .Call("_CVXR_LinOp__set_sparse_data", ptr, value, PACKAGE = "CVXR")
  get_sparse_data <- function() .Call("_CVXR_LinOp__get_sparse_data", ptr, PACKAGE = "CVXR")
  set_dense_data <- function(value) .Call("_CVXR_LinOp__set_dense_data", ptr, value, PACKAGE = "CVXR")
  get_dense_data <- function() .Call("_CVXR_LinOp__get_dense_data", ptr, PACKAGE = "CVXR")
  set_type <- function(value) {
    ##value <- match.arg(value, operatorType)
    ## Make zero based index!
    index <- match(value, operatorType) - 1L
    .Call("_CVXR_LinOp__set_type", ptr, index, PACKAGE = "CVXR")
  }
  get_type <- function() {
    index <- .Call("_CVXR_LinOp__get_type", ptr, PACKAGE = "CVXR")
    operatorType[index + 1L]
  }
  set_size <- function(value) .Call("_CVXR_LinOp__set_size", ptr, value, PACKAGE = "CVXR")
  get_size <- function() .Call("_CVXR_LinOp__get_size", ptr, PACKAGE = "CVXR")
  set_slice <- function(value) .Call("_CVXR_LinOp__set_slice", ptr, value, PACKAGE = "CVXR")
  get_slice <- function() .Call("_CVXR_LinOp__get_slice", ptr, PACKAGE = "CVXR")
  args_push_back <- function(linOp) .Call("_CVXR_LinOp__args_push_back", ptr, linOp$getXPtr(), PACKAGE = "CVXR")
  slice_push_back <- function(intVector) .Call("_CVXR_LinOp__slice_push_back", ptr, intVector, PACKAGE = "CVXR")
  getArgs <- function() args
  get_id <- function() .Call("_CVXR_LinOp__get_id", ptr, PACKAGE = "CVXR")
  size_push_back <- function(value) .Call("_CVXR_LinOp__size_push_back", ptr, value, PACKAGE = "CVXR")
  toString <- function() {
    sparse <- get_sparse()
    if (sparse) {
      data <- paste(get_sparse_data(), collapse=", ")
    } else {
      data <- paste(get_dense_data(), collapse=", ")
    }
    sprintf("LinOp(id=%s, type=%s, size=[%s], args=%s, sparse=%s, data=[%s])",
            get_id(),
            get_type(),
            paste(get_size(), collapse=", "),
            args$toString(),
            sparse,
            data)
  }
  print <- function() print(toString())

  list(getXPtr = getXPtr,
  set_sparse = set_sparse,
  get_sparse = get_sparse,
  set_sparse_data = set_sparse_data,
  get_sparse_data = get_sparse_data,
  set_dense_data = set_dense_data,
  get_dense_data = get_dense_data,
  set_type = set_type,
  get_type = get_type,
  set_size = set_size,
  get_size = get_size,
  set_slice = set_slice,
  get_slice = get_slice,
  args_push_back = args_push_back,
  slice_push_back = slice_push_back,
  getArgs = getArgs,
  get_id = get_id,
  size_push_back = size_push_back,
  toString = toString,
  print = print)
}
