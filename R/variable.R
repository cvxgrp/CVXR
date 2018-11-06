#'
#' The Variable class.
#'
#' This class represents an optimization variable.
#'
#' @slot shape The dimensions of the variable.
#' @slot name (Optional) A character string representing the name of the variable.
#' @slot var_id (Internal) A unique identification number used internally.
#' @name Variable-class
#' @aliases Variable
#' @rdname Variable-class
.Variable <- setClass("Variable", representation(shape = "NumORNULL", name = "character", var_id = "integer", value = "ConstVal"),
                                  prototype(shape = NULL, name = NA_character_, var_id = NA_integer_, value = NA_real_), 
                                  validity = function(object) {
                                    if(!is.na(object@value))
                                      stop("[Variable: validation] value is an internal slot and should not be set by the user")
                                    return(TRUE)
                                  }, contains = "Leaf")

#' @param shape The dimensions of the variable.
#' @param name (Optional) A character string representing the name of the variable.
#' @rdname Variable-class
#' @examples
#' x <- Variable(3, name = "x0") ## 3-int variable
#' y <- Variable(3, 3, name = "y0") # Matrix variable
#' as.character(y)
#' id(y)
#' is_positive(x)
#' is_negative(x)
#' size(y)
#' name(y)
#' value(y) <- matrix(1:9, nrow = 3)
#' value(y)
#' grad(y)
#' variables(y)
#' canonicalize(y)
#' @export
Variable <- function(shape = NULL, name = NA_character_) { .Variable(shape = shape, name = name) }

setMethod("initialize", "Variable", function(.Object, ..., shape = NULL, name = NA_character_, var_id = get_id(), value = NA_real_) {
  .Object@var_id <- var_id
  if(is.na(name))
    .Object@name <- sprintf("%s%s", VAR_PREFIX, .Object@id)
  else
    .Object@name <- name
  .Object@value <- NA_real_
  callNextMethod(.Object, ..., shape = shape)
})

setMethod("show", "Variable", function(object) {
  paste("Variable(", paste(shape(x), collapse = ", "), ")", sep = "")
})

#' @param x,object A \linkS4class{Variable} object.
#' @rdname Variable-class
setMethod("as.character", "Variable", function(x) {
  paste("Variable(", paste(shape(x), collapse = ", "), ")", sep = "")
})

#' @describeIn Variable The unique ID of the variable.
setMethod("id", "Variable", function(object) { object@id })

#' @describeIn Variable The name of the variable.
#' @export
setMethod("name", "Variable", function(object) { object@name })

# Set the value of the primal variable.
setMethod("save_value", "Variable", function(object, value) {
  value <- validate_val(object, value)
  object@primal_value <- value
  object
})

#' @describeIn Variable The value of the variable.
setMethod("value", "Variable", function(object) { object@primal_value })

#' @param value The value to assign to the primal variable.
#' @describeIn Variable Set the value of the primal variable.
setReplaceMethod("value", "Variable", function(object, value) {
  object <- save_value(object, value)
  object
})

#' @describeIn Variable The sub/super-gradient of the variable represented as a sparse matrix.
setMethod("grad", "Variable", function(object) {
  # TODO: Do not assume shape is 2-D.
  len <- size(object)
  result <- list(sparseMatrix(i = 1:len, j = 1:len, x = rep(1, len)))
  names(result) <- as.character(id(object))
  result
})

#' @describeIn Variable Returns itself as a variable.
setMethod("variables", "Variable", function(object) { list(object) })

#' @describeIn Variable The canonical form of the variable.
setMethod("canonicalize", "Variable", function(object) {
  obj <- create_var(shape(object), id(object))
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
#' @aliases Bool
#' @rdname Bool-class
.Bool <- setClass("Bool", contains = "Variable")

#' @param rows The number of rows in the variable.
#' @param cols The number of columns in the variable.
#' @param name (Optional) A character string representing the name of the variable.
#' @rdname Bool-class
#' @examples
#' x <- Bool(3, name = "indicator") ## Boolean 3-vector
#' y <- Bool(3, 3) ## Matrix boolean
#' name(x)
#' as.character(x)
#' canonicalize(y)
#' is_positive(x)
#' is_negative(y)
#' @export
Bool <- function(rows = 1, cols = 1, name = NA_character_) { .Bool(rows = rows, cols = cols, name = name) }

setMethod("show", "Bool", function(object) {
  size <- size(object)
  cat("Bool(", size[1], ", ", size[2], ")", sep = "")
})

#' @param x,object A \linkS4class{Bool} object.
#' @rdname Bool-class
setMethod("as.character", "Bool", function(x) {
  size <- size(x)
  paste("Bool(", size[1], ", ", size[2], ")", sep = "")
})

#' @describeIn Bool Enforce that the variable be boolean.
setMethod("canonicalize", "Bool", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constr <- canon[[2]]
  list(obj, c(constr, list(BoolConstr(obj))))
})

#' @describeIn Bool A boolean variable is always positive or zero.
setMethod("is_positive", "Bool", function(object) { TRUE })

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
#' @aliases Int
#' @rdname Int-class
.Int <- setClass("Int", contains = "Variable")

#' @param rows The number of rows in the variable.
#' @param cols The number of columns in the variable.
#' @param name (Optional) A character string representing the name of the variable.
#' @examples
#' x <- Int(3, name = "i") ## 3-int variable
#' y <- Int(3, 3, name = "j") # Matrix variable
#' as.character(y)
#' id(y)
#' is_positive(x)
#' is_negative(x)
#' size(y)
#' name(y)
#' value(y) <- matrix(1:9, nrow = 3)
#' value(y)
#' grad(y)
#' variables(y)
#' canonicalize(y)
#' @rdname Int-class
#' @export
Int <- function(rows = 1, cols = 1, name = NA_character_) { .Int(rows = rows, cols = cols, name = name) }

setMethod("show", "Int", function(object) {
  size <- size(object)
  cat("Int(", size[1], ", ", size[2], ")", sep = "")
})

#' @param x,object An \linkS4class{Int} object.
#' @rdname Int-class
setMethod("as.character", "Int", function(x) {
  size <- size(x)
  paste("Int(", size[1], ", ", size[2], ")", sep = "")
})

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
#' @aliases NonNegative
#' @rdname NonNegative-class
.NonNegative <- setClass("NonNegative", contains = "Variable")

#' @param rows The number of rows in the variable.
#' @param cols The number of columns in the variable.
#' @param name (Optional) A character string representing the name of the variable.
#' @rdname NonNegative-class
#' @examples
#' x <- NonNegative(3, 3)
#' as.character(x)
#' canonicalize(x)
#' is_positive(x)
#' is_negative(x)
#' @export
NonNegative <- function(rows = 1, cols = 1, name = NA_character_) { .NonNegative(rows = rows, cols = cols, name = name) }

setMethod("show", "NonNegative", function(object) {
  size <- size(object)
  cat("NonNegative(", size[1], ", ", size[2], ")", sep = "")
})

#' @param x,object A \linkS4class{NonNegative} object.
#' @rdname NonNegative-class
setMethod("as.character", "NonNegative", function(x) {
  size <- size(x)
  paste("NonNegative(", size[1], ", ", size[2], ")", sep = "")
})

#' @describeIn NonNegative Enforce that the variable be non-negative.
setMethod("canonicalize", "NonNegative", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constr <- canon[[2]]
  list(obj, c(constr, list(create_geq(obj))))
})

#' @describeIn NonNegative Always true since the variable is non-negative.
setMethod("is_positive", "NonNegative", function(object) { TRUE })

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
#' @aliases SemidefUpperTri
#' @rdname SemidefUpperTri-class
.SemidefUpperTri <- setClass("SemidefUpperTri", representation(n = "numeric"), contains = "Variable")

#' @param n The number of rows/columns in the matrix.
#' @param name (Optional) A character string representing the name of the variable.
#' @rdname SemidefUpperTri-class
#' @examples
#' x <- SemidefUpperTri(3)
#' as.character(x)
#' get_data(x)
#' canonicalize(x)
#' @export
SemidefUpperTri <- function(n, name = NA_character_) { .SemidefUpperTri(n = n, name = name) }

setMethod("initialize", "SemidefUpperTri", function(.Object, ..., rows, cols, name = NA_character_, n) {
  .Object@n = n
  callNextMethod(.Object, ..., rows = n*(n+1)/2, cols = 1, name = name)
})

setMethod("show", "SemidefUpperTri", function(object) {
  cat("SemidefUpperTri(", object@n, ")", sep = "")
})

#' @param x,object A \linkS4class{SemidefUpperTri} object.
#' @rdname SemidefUpperTri-class
setMethod("as.character", "SemidefUpperTri", function(x) {
  paste("SemidefUpperTri(", x@n, ")", sep = "")
})

#' @describeIn SemidefUpperTri Returns \code{list(n, name)}.
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
#' @examples
#' x <- Semidef(5) ## 5 by 5 semidefinite matrix expression
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
#' @aliases SymmetricUpperTri
#' @rdname SymmetricUpperTri-class
.SymmetricUpperTri <- setClass("SymmetricUpperTri", representation(n = "numeric"), contains = "Variable")

#' @param n The number of rows/columns in the matrix.
#' @param name (Optional) A character string representing the name of the variable.
#' @rdname SymmetricUpperTri-class
#' @examples
#' x <- SymmetricUpperTri(3, name="s3")
#' name(x)
#' get_data(x)
#' @export
SymmetricUpperTri <- function(n, name = NA_character_) { .SymmetricUpperTri(n = n, name = name) }

setMethod("initialize", "SymmetricUpperTri", function(.Object, ..., rows, cols, name = NA_character_, n) {
  .Object@n = n
  callNextMethod(.Object, ..., rows = n*(n+1)/2, cols = 1, name = name)
})

setMethod("show", "SymmetricUpperTri", function(object) {
  cat("SymmetricUpperTri(", object@n, ")", sep = "")
})

#' @param x,object A \linkS4class{SymmetricUpperTri} object.
#' @rdname SymmetricUpperTri-class
setMethod("as.character", "SymmetricUpperTri", function(x) {
  paste("SymmetricUpperTri(", x@n, ")", sep = "")
})

#' @describeIn SymmetricUpperTri Returns \code{list(n, name)}.
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
#' @examples
#' x <- Symmetric(3, name="s3")
#' @export
Symmetric <- function(n, name = NA_character_) {
  var <- SymmetricUpperTri(n, name)
  fill_mat <- Constant(upper_tri_to_full(n))
  Reshape(fill_mat %*% var, as.integer(n), as.integer(n))
}
