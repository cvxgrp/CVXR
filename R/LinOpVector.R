#' Simple object system that backs the C++ LinOpVector class
#' We need this object because we don't want items pushed on to a C++ vector disposed of!
#' So we need to hold on to a handle for pushed objects here.
#' @param ... named args to object 
make_LinOpVector <- function(...) {
  self <- environment()
  class(self) <- c("LinOpVector", class(self))
  self$self <- self  

  ## Create C++ object
  self$xptr <- .Call(paste0("_CVXR_LinOpVector__new"), PACKAGE = "CVXR")
  self$lin_ops <- list()
  
  ## Define Methods
  self$push_back <- function(lin_op) { # lin_op is the R object backing the C++ object.
    n <- length(self$lin_ops)
    self$lin_ops[[n + 1]] <- lin_op
    ## Needs modification by hand for arguments
    .Call("_CVXR_LinOpVector__push_back", self$xptr , lin_op$xptr, PACKAGE = "CVXR")
  }
  self
}

#' @method print LinOpVector
print.LinOpVector <- function(x, ...) {
  comps <- ls(x)
  methods <- sapply(comps, function(name) is.function(x[[name]]))
  others <- comps[!methods]
  out <- sprintf("LinOpVector[%d] class", length(x$lin_ops))
  if (length(others) > 0) {
    out <- c(
      out,
      " " = paste0("Slots: ", paste(others, collapse = ", "))
    )
  }
  if (length(methods) > 0) {
    out <- c(
      out,
      " " = paste0("Methods: ", paste(comps[methods], collapse = ", "))
    )
  }
  cli::cli_bullets(out)
}

#' This differs from LinOpVector above in only that the vector elements are const *, so that they cannot be modified.
#' 
#' @param ... named args to object 
make_ConstLinOpVector <- function(...) {
  self <- environment()
  class(self) <- c("ConstLinOpVector", class(self))
  self$self <- self  

  ## Create C++ object
  self$xptr <- .Call(paste0("_CVXR_LinOpVector__new"), PACKAGE = "CVXR")
  self$lin_ops <- list()
  
  ## Define Methods
  self$push_back <- function(lin_op) { # lin_op is the R object backing the C++ object.
    n <- length(self$lin_ops)
    self$lin_ops[[n + 1]] <- lin_op
    ## Needs modification by hand for arguments
    .Call("_CVXR_LinOpVector__push_back", self$xptr , lin_op$xptr, PACKAGE = "CVXR")
  }
  self
}

#' @method print ConstLinOpVector
print.ConstLinOpVector <- function(x, ...) {
  comps <- ls(x)
  methods <- sapply(comps, function(name) is.function(x[[name]]))
  others <- comps[!methods]
  out <- sprintf("ConstLinOpVector[%d] class", length(x$lin_ops))
  if (length(others) > 0) {
    out <- c(
      out,
      " " = paste0("Slots: ", paste(others, collapse = ", "))
    )
  }
  if (length(methods) > 0) {
    out <- c(
      out,
      " " = paste0("Methods: ", paste(comps[methods], collapse = ", "))
    )
  }
  cli::cli_bullets(out)
}

