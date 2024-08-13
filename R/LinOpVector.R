## LinOpVector class shadowing CPP class
CVXcanon.LinOpVector <- function() {
  linOps <- list()
  ptr <- .Call("_CVXR_LinOpVector__new", PACKAGE = "CVXR")
  getXPtr <- function() ptr
  getList <- function() linOps
  push_back <- function(linOp) {
    n <- length(linOps)
    linOps[[n + 1L]] <- linOp
    ## Needs modification by hand for arguments
    .Call("_CVXR_LinOpVector__push_back", ptr , linOp$getXPtr(), PACKAGE = "CVXR")
  }
  toString <- function() {
    result <- sapply(linOps, function(x) x$toString())
    result <- paste(result, collapse = ", ")
    sprintf("[ %s ]", result)
  }
  print <- function() print(self$toString())
  list(getXPtr = getXPtr,
       getList = getList,
       push_back = push_back,
       toString = toString,
       print = print)
}
