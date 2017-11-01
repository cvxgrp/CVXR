#'
#' The Variable class.
#'
#' This class represents an optimization variable.
#' 
#' @slot id (Internal) A unique identification number used internally.
#' @slot rows The number of rows in the variable.
#' @slot cols The number of columns in the variable.
#' @slot name (Optional) A character string representing the name of the variable.
#' @slot primal_value (Internal) The primal value of the variable stored internally.
#' @name Variable-class
#' @rdname Variable-class
#' @export
.Variable <- setClass("Variable", representation(id = "integer", rows = "numeric", cols = "numeric", name = "character", primal_value = "ConstVal"),
                                 prototype(rows = 1, cols = 1, name = NA_character_, primal_value = NA_real_),
                                 validity = function(object) {
                                   if(!is.na(object@primal_value))
                                     stop("[Variable: validation] primal_value is an internal slot and should not be set by user")
                                   return(TRUE)
                                 }, contains = "Leaf")

#'
#' Variable Constructor
#'
#' Construct a \linkS4class{Variable} class object.
#' 
#' @param rows The number of rows in the variable.
#' @param cols The number of columns in the variable.
#' @param name (Optional) A character string representing the name of the variable.
#' @return A \linkS4class{Variable} object.
#' @name Variable
#' @rdname Variable-class
#' @export
Variable <- function(rows = 1, cols = 1, name = NA_character_) { .Variable(rows = rows, cols = cols, name = name) }

#'
#' Variable Initialization
#'
#' @name Variable
#' @rdname Variable-class
setMethod("initialize", "Variable", function(.Object, ..., id = get_id(), rows = 1, cols = 1, name = NA_character_, primal_value = NA_real_) {
  .Object@id <- id
  .Object@rows <- rows
  .Object@cols <- cols
  if(is.na(name))
    .Object@name <- sprintf("%s%s", VAR_PREFIX, .Object@id)
  else
    .Object@name <- name
  .Object@primal_value <- primal_value
  callNextMethod(.Object, ...)
})

setMethod("show", "Variable", function(object) {
  size <- size(object)
  cat("Variable(", size[1], ", ", size[2], ")", sep = "")
})

setMethod("as.character", "Variable", function(x) {
  size <- size(x)
  paste("Variable(", size[1], ", ", size[2], ")", sep = "")
})

#' @describeIn Variable The unique ID of the variable.
setMethod("id", "Variable", function(object) { object@id })

#' @rdname sign-methods
#' @describeIn Variable A logical value indicating whether the variable is positive.
setMethod("is_positive", "Variable", function(object) { FALSE })

#' @rdname sign-methods
#' @describeIn Variable A logical value indicating whether the variable is negative.
setMethod("is_negative", "Variable", function(object) { FALSE })

#' @rdname size
#' @describeIn Variable The \code{c(row, col)} dimensions of the variable.
setMethod("size", "Variable", function(object) { c(object@rows, object@cols) })

setMethod("get_data", "Variable", function(object) { list(object@rows, object@cols, object@name) })

#' @rdname name
#' @describeIn Variable The name of the variable.
setMethod("name", "Variable", function(object) { object@name })

# Set the value of the primal variable.
setMethod("save_value", "Variable", function(object, value) {
  value <- validate_val(object, value)
  object@primal_value <- value
  object
})

#' @rdname value-methods
#' @describeIn Variable The value of the variable.
setMethod("value", "Variable", function(object) { object@primal_value })

#' @rdname value-methods
#' @describeIn Variable Set the value of the primal variable.
setReplaceMethod("value", "Variable", function(object, value) {
  object <- save_value(object, value)
  object
})

#' @rdname grad
#' @describeIn Variable The sub/super-gradient of the variable represented as a sparse matrix.
setMethod("grad", "Variable", function(object) {
  len <- prod(size(object))
  result <- list(sparseMatrix(i = 1:len, j = 1:len, x = rep(1, len)))
  names(result) <- as.character(object@id)
  result
})

#' @rdname expression-parts
#' @describeIn Variable Returns itself as a variable.
setMethod("variables", "Variable", function(object) { list(object) })

#' @rdname canonicalize
#' @describeIn Variable The canonical form of the variable.
setMethod("canonicalize", "Variable", function(object) {
  obj <- create_var(size(object), object@id)
  list(obj, list())
})

#'
#' The Bool class.
#' 
#' This class represents a boolean variable.
#'
#' @slot id (Internal) A unique identification number used internally.
#' @slot rows The number of rows in the variable.
#' @slot cols The number of columns in the variable.
#' @slot name (Optional) A character string representing the name of the variable.
#' @slot primal_value (Internal) The primal value of the variable stored internally.
#' @name Bool-class
#' @rdname Bool-class
#' @export
.Bool <- setClass("Bool", contains = "Variable")

#'
#' Boolean Variable Constructor
#'
#' Construct a \linkS4class{Bool} class object.
#' 
#' @param rows The number of rows in the variable.
#' @param cols The number of columns in the variable.
#' @param name (Optional) A character string representing the name of the variable.
#' @return A \linkS4class{Bool} object.
#' @name Bool
#' @rdname Bool-class
#' @export
Bool <- function(rows = 1, cols = 1, name = NA_character_) { .Bool(rows = rows, cols = cols, name = name) }

setMethod("show", "Bool", function(object) {
  size <- size(object)
  cat("Bool(", size[1], ", ", size[2], ")", sep = "")
})

setMethod("as.character", "Bool", function(x) {
  size <- size(x)
  paste("Bool(", size[1], ", ", size[2], ")", sep = "")
})

#' @rdname canonicalize
#' @describeIn Bool Enforce that the variable be boolean.
setMethod("canonicalize", "Bool", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constr <- canon[[2]]
  list(obj, c(constr, list(BoolConstr(obj))))
})

#' @rdname sign-methods
#' @describeIn Bool A boolean variable is always positive or zero.
setMethod("is_positive", "Bool", function(object) { TRUE })

#' @rdname sign-methods
#' @describeIn Bool A boolean variable is never negative.
setMethod("is_negative", "Bool", function(object) { FALSE })

#'
#' The Int class.
#' 
#' This class represents an integer variable.
#'
#' @slot id (Internal) A unique identification number used internally.
#' @slot rows The number of rows in the variable.
#' @slot cols The number of columns in the variable.
#' @slot name (Optional) A character string representing the name of the variable.
#' @slot primal_value (Internal) The primal value of the variable stored internally.
#' @name Int-class
#' @rdname Int-class
#' @export
.Int <- setClass("Int", contains = "Variable")

#'
#' Integer Variable Constructor
#'
#' Construct a \linkS4class{Int} class object.
#' 
#' @param rows The number of rows in the variable.
#' @param cols The number of columns in the variable.
#' @param name (Optional) A character string representing the name of the variable.
#' @return A \linkS4class{Int} object.
#' @name Int
#' @rdname Int-class
#' @export
Int <- function(rows = 1, cols = 1, name = NA_character_) { .Int(rows = rows, cols = cols, name = name) }

setMethod("show", "Int", function(object) {
  size <- size(object)
  cat("Int(", size[1], ", ", size[2], ")", sep = "")
})

setMethod("as.character", "Int", function(x) {
  size <- size(x)
  paste("Int(", size[1], ", ", size[2], ")", sep = "")
})

#' @rdname canonicalize
#' @describeIn Int Enforce that the variable be an integer.
setMethod("canonicalize", "Int", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constr <- canon[[2]]
  list(obj, c(constr, list(IntConstr(obj))))
})

#'
#' The NonNegative class.
#'
#' This class represents a variable constrained to be non-negative.
#' 
#' @slot id (Internal) A unique identification number used internally.
#' @slot rows The number of rows in the variable.
#' @slot cols The number of columns in the variable.
#' @slot name (Optional) A character string representing the name of the variable.
#' @slot primal_value (Internal) The primal value of the variable stored internally.
#' @name NonNegative-class
#' @rdname NonNegative-class
#' @export
.NonNegative <- setClass("NonNegative", contains = "Variable")

#'
#' Non-Negative Variable Constructor
#'
#' Construct a \linkS4class{NonNegative} class object.
#' 
#' @param rows The number of rows in the variable.
#' @param cols The number of columns in the variable.
#' @param name (Optional) A character string representing the name of the variable.
#' @return A \linkS4class{NonNegative} object.
#' @name NonNegative
#' @rdname NonNegative-class
#' @export
NonNegative <- function(rows = 1, cols = 1, name = NA_character_) { .NonNegative(rows = rows, cols = cols, name = name) }

setMethod("show", "NonNegative", function(object) {
  size <- size(object)
  cat("NonNegative(", size[1], ", ", size[2], ")", sep = "")
})

setMethod("as.character", "NonNegative", function(x) {
  size <- size(x)
  paste("NonNegative(", size[1], ", ", size[2], ")", sep = "")
})

#' @rdname canonicalize
#' @describeIn NonNegative Enforce that the variable be non-negative.
setMethod("canonicalize", "NonNegative", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constr <- canon[[2]]
  list(obj, c(constr, list(create_geq(obj))))
})

#' @rdname sign-methods
#' @describeIn NonNegative Always true since the variable is non-negative.
setMethod("is_positive", "NonNegative", function(object) { TRUE })

#' @rdname sign-methods
#' @describeIn NonNegative Always false since the variable is non-negative.
setMethod("is_negative", "NonNegative", function(object) { FALSE })

#'
#' The SemidefUpperTri class.
#'
#' This class represents the upper triangular part of a positive semidefinite variable.
#'
#' @slot id (Internal) A unique identification number used internally.
#' @slot n The number of rows/columns in the matrix.
#' @slot rows The number of rows in the variable.
#' @slot cols The number of columns in the variable.
#' @slot name (Optional) A character string representing the name of the variable.
#' @slot primal_value (Internal) The primal value of the variable stored internally.
#' @name SemidefUpperTri-class
#' @rdname SemidefUpperTri-class
.SemidefUpperTri <- setClass("SemidefUpperTri", representation(n = "numeric"), contains = "Variable")

#'
#' Upper Triangle of Positive Semidefinite Variable
#' 
#' Construct a \linkS4class{SemidefUpperTri} class object.
#' 
#' @param n The number of rows/columns in the matrix.
#' @param name (Optional) A character string representing the name of the variable.
#' @return A \linkS4class{SemidefUpperTri} object.
#' @name SemidefUpperTri
#' @rdname SemidefUpperTri-class
SemidefUpperTri <- function(n, name = NA_character_) { .SemidefUpperTri(n = n, name = name) }

#' @name SemidefUpperTri
#' @rdname SemidefUpperTri-class
setMethod("initialize", "SemidefUpperTri", function(.Object, ..., rows, cols, name = NA_character_, n) {
  .Object@n = n
  callNextMethod(.Object, ..., rows = n*(n+1)/2, cols = 1, name = name)
})

setMethod("show", "SemidefUpperTri", function(object) {
  cat("SemidefUpperTri(", object@n, ")", sep = "")
})

setMethod("as.character", "SemidefUpperTri", function(x) {
  paste("SemidefUpperTri(", x@n, ")", sep = "")
})

setMethod("get_data", "SemidefUpperTri", function(object) { list(object@n, object@name) })

#
# Upper Triangle to Full Matrix
#
# Returns a coefficient matrix to create a symmetric matrix.
#
# @param n The width/height of the matrix
# @return The coefficient matrix.
upper_tri_to_full <- function(n) {
  if(n == 0)
    return(sparseMatrix(i = c(), j = c(), dims = c(0, 0)))
  
  entries <- floor(n*(n+1)/2)
  
  val_arr <- c()
  row_arr <- c()
  col_arr <- c()
  count <- 1
  for(i in 1:n) {
    for(j in i:n) {
      # Index in the original matrix
      col_arr <- c(col_arr, count)
      
      # Index in the filled matrix
      row_arr <- c(row_arr, (j-1)*n + i)
      val_arr <- c(val_arr, 1.0)
      if(i != j) {
        # Index in the original matrix
        col_arr <- c(col_arr, count)
        
        # Index in the filled matrix
        row_arr <- c(row_arr, (i-1)*n + j)
        val_arr <- c(val_arr, 1.0)
      }
      count <- count + 1
    }
  }
  sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = c(n^2, entries))
}

#' @describeIn SemidefUpperTri Enforce that the variable be positive semidefinite.
setMethod("canonicalize", "SemidefUpperTri", function(object) {
  # Variable must be semidefinite and symmetric
  upper_tri <- create_var(c(size(object)[1], 1), object@id)
  fill_coef <- upper_tri_to_full(object@n)
  fill_coef <- create_const(fill_coef, c(object@n^2, size(object)[1]), sparse = TRUE)
  full_mat = lo.mul_expr(fill_coef, upper_tri, c(object@n^2, 1))
  full_mat <- lo.reshape(full_mat, c(object@n, object@n))
  list(upper_tri, list(SDP(full_mat, enforce_sym = FALSE)))
})

#'
#' Positive Semidefinite Variable
#'
#' An expression representing a positive semidefinite matrix.
#'
#' @param n The number of rows/columns in the matrix.
#' @param name (Optional) A character string representing the name of the variable.
#' @return An \linkS4class{Expression} representing the positive semidefinite matrix.
#' @rdname Semidef
#' @export
Semidef <- function(n, name = NA_character_) {
  var <- SemidefUpperTri(n, name)
  fill_mat <- Constant(upper_tri_to_full(n))
  Reshape(fill_mat %*% var, n, n)
}

#'
#' The SymmetricUpperTri class.
#'
#' This class represents the upper triangular part of a symmetric variable.
#'
#' @slot id (Internal) A unique identification number used internally.
#' @slot n The number of rows/columns in the matrix.
#' @slot rows The number of rows in the variable.
#' @slot cols The number of columns in the variable.
#' @slot name (Optional) A character string representing the name of the variable.
#' @slot primal_value (Internal) The primal value of the variable stored internally.
#' @name SymmetricUpperTri-class
#' @rdname SymmetricUpperTri-class
.SymmetricUpperTri <- setClass("SymmetricUpperTri", representation(n = "numeric"), contains = "Variable")

#'
#' Upper Triangle of Symmetric Variable
#' 
#' Construct a \linkS4class{SymmetricUpperTri} class object.
#' 
#' @param n The number of rows/columns in the matrix.
#' @param name (Optional) A character string representing the name of the variable.
#' @return A \linkS4class{SymmetricUpperTri} object.
#' @name SymmetricUpperTri
#' @rdname SymmetricUpperTri-class
SymmetricUpperTri <- function(n, name = NA_character_) { .SymmetricUpperTri(n = n, name = name) }

#' @name SymmetricUpperTri
#' @rdname SymmetricUpperTri-class
setMethod("initialize", "SymmetricUpperTri", function(.Object, ..., rows, cols, name = NA_character_, n) {
  .Object@n = n
  callNextMethod(.Object, ..., rows = n*(n+1)/2, cols = 1, name = name)
})

setMethod("show", "SymmetricUpperTri", function(object) {
  cat("SymmetricUpperTri(", object@n, ")", sep = "")
})

setMethod("as.character", "SymmetricUpperTri", function(x) {
  paste("SymmetricUpperTri(", x@n, ")", sep = "")
})

setMethod("get_data", "SymmetricUpperTri", function(object) { list(object@n, object@name) })

#' @describeIn SemidefUpperTri Enforce that the variable be symmetric.
setMethod("canonicalize", "SymmetricUpperTri", function(object) {
  upper_tri <- create_var(c(size(object)[1], 1), object@id)
  list(upper_tri, list())
})

#'
#' Symmetric Variable
#'
#' An expression representing a symmetric matrix.
#'
#' @param n The number of rows/columns in the matrix.
#' @param name (Optional) A character string representing the name of the variable.
#' @return An \linkS4class{Expression} representing the symmetric matrix.
#' @rdname Symmetric
#' @export
Symmetric <- function(n, name = NA_character_) {
  var <- SymmetricUpperTri(n, name)
  fill_mat <- Constant(upper_tri_to_full(n))
  Reshape(fill_mat %*% var, as.integer(n), as.integer(n))
}
