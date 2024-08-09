## LinOp class shadowing CPP class
cvxcore.LinOp <- R6::R6Class(
  classname = "cvxcore.LinOp",  ## We don't use the class name in any meaningful way and NULL speeds up things
  private = list(
    args = NA,
    ptr = NA
  ),
  active = list(
    sparse = function(value) {
      if (missing(value)) {
        .Call("_CVXR_LinOp__get_sparse", private$ptr, PACKAGE = "CVXR")
      } else {
        ## value should be a boolean
        .Call("_CVXR_LinOp__set_sparse",  private$ptr, value, PACKAGE = "CVXR")
      }
    }
   ,
    sparse_data = function(value) {
      if (missing(value)) {
        .Call("_CVXR_LinOp__get_sparse_data", private$ptr, PACKAGE = "CVXR")
      } else {
        ## value should be a dgCMatrix-class
        .Call("_CVXR_LinOp__set_sparse_data", private$ptr, value, PACKAGE = "CVXR")
      }
    }
   ,
    dense_data = function(value) {
      if (missing(value)) {
        .Call("_CVXR_LinOp__get_dense_data", private$ptr, PACKAGE = "CVXR")
      } else {
        ## value should be a matrix
        .Call("_CVXR_LinOp__set_dense_data", private$ptr, value, PACKAGE = "CVXR")
      }
    }
   ,
    type = function(value) {
      if (missing(value)) {
        index <- .Call("_CVXR_LinOp__get_type", private$ptr, PACKAGE = "CVXR")
        ## make 1-based index
        LINOP_TYPES[index]
      } else {
        index <- match(value, LINOP_TYPES)
        .Call("_CVXR_LinOp__set_type", private$ptr, index, PACKAGE = "CVXR")
      }
    }
   ,
    data_ndim = function(value) {
      if (missing(value)) {
        .Call("_CVXR_LinOp__get_data_ndim", private$ptr, PACKAGE = "CVXR")
      } else {
        ## value is an integer vector
        .Call("_CVXR_LinOp__set_data_ndim", private$ptr, value, PACKAGE = "CVXR")
      }
    }
  ),
  public = list(
    initialize = function(type, size, args) {
      private$args = R6List$new()
      ## Create a new LinOp on the C side
      private$ptr <- .Call("_CVXR_LinOp__new", type, dim, args, PACKAGE = "CVXR")
      ## Associate args on R side with the args on the C side.
      if (!is.null(type)) {
        self$type <- type
      }
      if (!is.null(size)) {
        self$size <- size
      }
      if (!is.null(args)) {
        for (x in args) self$args_push_back(x)
      }
      if (!is.null(data)) {
        self$dense_data <- data
      }
    }
   ,
    args_push_back = function(R6LinOp) {
      private$args$append(R6LinOp)
      .Call("_CVXR_LinOp__args_push_back", private$ptr, R6LinOp$getXPtr(), PACKAGE = "CVXR")
    }
   ,
    slice_push_back = function(anIntVector) {
      .Call("_CVXR_LinOp__slice_push_back", private$ptr,
            anIntVector, PACKAGE = "CVXR")
    }
   ,
    getXPtr = function() {
      private$ptr
    }
   ,
    getArgs = function() {
      private$args
    }
   ,
    get_id = function() {
      .Call("_CVXR_LinOp__get_id", private$ptr, PACKAGE = "CVXR")
    }
   ,
    size_push_back = function(value) {
      .Call("_CVXR_LinOp__size_push_back", private$ptr, value, PACKAGE = "CVXR")
    }
   ,
    toString = function() {
      sparse <- self$sparse
      if (sparse) {
        data <- paste(self$sparse_data, collapse=", ")
      } else {
        data <- paste(self$dense_data, collapse=", ")
      }
      sprintf("LinOp(id=%s, type=%s, size=[%s], args=%s, sparse=%s, data=[%s])",
              self$get_id(),
              self$type,
              paste(self$size, collapse=", "),
              private$args$toString(),
              sparse,
              data)
    }
   ,
    print = function() {
      print(self$toString())
    }
    
  ))
