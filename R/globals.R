#' @importFrom utils assignInMyNamespace
#' @importFrom R.utils intToBin
#' @importFrom Rmpfr mpfr getPrec
#' @importFrom gmp as.bigq as.bigz is.bigq is.bigz is.whole numerator denominator asNumeric
#' @importFrom bit64 as.integer64 as.bitstring
#' @importClassesFrom gmp bigq bigz

##
## Global variable for package CVXR
##
.CVXR.options <- list(idCounter = 0L)

#'
#' Set ID Counter
#'
#' Set the CVXR variable/constraint identification number counter.
#'
#' @param value The value to assign as ID.
#' @return the changed value of the package global \code{.CVXR.options}.
#' @export
#' @examples
#' \dontrun{
#'   setIdCounter(value = 0L)
#' }
setIdCounter <- function(value = 0L) {
    .CVXR.options$idCounter <- value
    assignInMyNamespace(".CVXR.options", .CVXR.options)
    .CVXR.options
}

#'
#' Reset Options
#'
#' Reset the global package variable \code{.CVXR.options}.
#'
#' @return The default value of CVXR package global \code{.CVXR.options}.
#' @export
#' @examples
#' \dontrun{
#'   resetOptions()
#' }
resetOptions <- function() {
    assignInMyNamespace(".CVXR.options", list(idCounter = 0L))
    .CVXR.options
}

#'
#' Get ID
#'
#' Get the next identifier value.
#'
#' @return A new unique integer identifier.
#' @export
#' @examples
#' \dontrun{
#'    get_id()
#' }
get_id <- function() {
    id <- .CVXR.options$idCounter <- .CVXR.options$idCounter + 1L
    assignInMyNamespace(".CVXR.options", .CVXR.options)
    id
}
