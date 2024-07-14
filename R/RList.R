#' Simple object system that will be sufficient to backing C++ classes for us.
#' @param ... named args to object 
make_RList <- function(...) {
  self <- environment()
  class(self) <- c("RList", class(self))
  self$self <- self
  self$rlist <- list()
  ## Methods
  self$append <- function(what) {
    n <- length(self$rlist)
    self$rlist[[n + 1]] <- what
  }
  self
}

#' @method print RList
print.RList <- function(x, ...) {
  comps <- ls(x)
  methods <- sapply(comps, function(name) is.function(x[[name]]))
  others <- comps[!methods]
  out <- sprintf("RListin[%d] class", length(x$lin_ops))
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
