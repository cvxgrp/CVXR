##
## Global variable for package cvxr
##
.cvxr.options <- list(idCounter = 0L)

#' Set the cvxr idCounter
#'
#' @param value, the value to assign
#' @return the changed value of the package global \code{.cvxr.options}
#' @export
#' @examples
#' \dontrun{
#'   setIdCounter(value = 0L)
#' }
#'
setIdCounter <- function(value = 0L) {
    .cvxr.options$idCounter <- value
    assignInMyNamespace(".cvxr.options", .cvxr.options)
    .cvxr.options
}

#' Reset the global package variable \code{.cvxr.options}
#'
#' @return the default value value of cvxr package global \code{.cvxr.options}
#' @export
#' @examples
#' \dontrun{
#'   resetOptions()
#' }
#'
resetOptions <- function() {
    assignInMyNamespace(".cvxr.options", list(idCounter = 0L))
    .cvxr.options
}

#' Get the next identifier value
#'
#' @return a new unique integer identifier
#' @export
#' @examples
#' \dontrun{
#'    get_id()
#' }
#'
get_id <- function() {
    id <- .cvxr.options$idCounter <- .cvxr.options$idCounter + 1L
    assignInMyNamespace(".cvxr.options", .cvxr.options)
    id
}
