#' @importFrom utils assignInMyNamespace
#' @importFrom R.utils intToBin
#' @importFrom Rmpfr mpfr getPrec
#' @importFrom gmp as.bigq as.bigz is.bigq is.bigz is.whole numerator denominator asNumeric
#' @importFrom bit64 as.integer64 as.bitstring
#' @importFrom reticulate import import_from_path py_module_available
#' @importClassesFrom gmp bigq bigz

##
## Global variable for package CVXR
## idCounter, numpy handle, scipy.sparse handle
##
.CVXR.options <- list(idCounter = 0L, np = NULL, sp = NULL, mosekglue = NULL)

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
    assignInMyNamespace(".CVXR.options", list(idCounter = 0L, np = NULL, sp = NULL, mosekglue = NULL))
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

#'
#' Get scipy handle
#'
#' Get the scipy handle or fail if not available
#'
#' @return the scipy handle
#' @export
#' @examples
#' \dontrun{
#'    get_sp
#' }
get_sp <- function() {
    sp <- .CVXR.options$sp
    if (is.null(sp)) {
        stop("Scipy not available")
    }
    sp
}

#'
#' Get numpy handle
#'
#' Get the numpy handle or fail if not available
#'
#' @return the numpy handle
#' @export
#' @examples
#' \dontrun{
#'    get_np
#' }
get_np <- function() {
    np <- .CVXR.options$np
    if (is.null(np)) {
        stop("Numpy not available")
    }
    np
}

#'
#' Get our mosekglue handle
#'
#' Get the mosekglue handle or fail if not available
#'
#' @return the mosekglue handle
#' @export
#' @examples
#' \dontrun{
#'    get_mosekglue
#' }
get_mosekglue <- function() {
    mosekglue <- .CVXR.options$mosekglue
    if (is.null(mosekglue)) {
        stop("CVXR python mosekglue not available")
    }
    mosekglue
}

#'
#' Get our gurobiglue handle
#'
#' Get the gurobiglue handle or fail if not available
#'
#' @return the gurobiglue handle
#' @export
#' @examples
#' \dontrun{
#'    get_gurobiglue
#' }
get_gurobiglue <- function() {
    gurobiglue <- .CVXR.options$gurobiglue
    if (is.null(gurobiglue)) {
        stop("CVXR python gurobiglue not available")
    }
    gurobiglue
}

