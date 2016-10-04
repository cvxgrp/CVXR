.Variable <- setClass("Variable", representation(id = "integer", rows = "numeric", cols = "numeric", name = "character", primal_value = "ConstVal"),
                                 prototype(id = get_id(), rows = 1, cols = 1, name = NA_character_, primal_value = NA_real_),
                                 validity = function(object) {
                                   if(!is.na(object@primal_value))
                                     stop("[Variable: validation] primal_value is an internal slot and should not be set by user")
                                   return(TRUE)
                                 }, contains = "Leaf")

Variable <- function(rows = 1, cols = 1, name = NA_character_) { .Variable(rows = rows, cols = cols, name = name) }

setMethod("initialize", "Variable", function(.Object, ..., id = get_id(), rows = 1, cols = 1, name = NA_character_, primal_value = NA_real_) {
  .Object@rows <- rows
  .Object@cols <- cols
  .Object@id <- id
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

setMethod("is_positive", "Variable", function(object) { FALSE })
setMethod("is_negative", "Variable", function(object) { FALSE })
setMethod("size", "Variable", function(object) { c(object@rows, object@cols) })
setMethod("get_data", "Variable", function(object) { list(object@rows, object@cols, object@name) })
setMethod("name", "Variable", function(object) { object@name })

setMethod("save_value", "Variable", function(object, value) {
  object@primal_value <- value
  object
})
setMethod("value", "Variable", function(object) { object@primal_value })
setReplaceMethod("value", "Variable", function(object, value) {
  value <- validate_val(object, value)
  object <- save_value(object, value)
  object
})
setMethod("grad", "Variable", function(object) {
  len <- prod(size(object))
  result <- list(sparseMatrix(i = 1:len, j = 1:len, x = rep(1, len)))
  names(result) <- as.character(object@id)
  result
})
setMethod("variables", "Variable", function(object) { list(object) })
setMethod("canonicalize", "Variable", function(object) {
  obj <- create_var(size(object), object@id)
  list(obj, list())
})

# Boolean variable
.Bool <- setClass("Bool", contains = "Variable")
Bool <- function(rows = 1, cols = 1, name = NA_character_) { .Bool(rows = rows, cols = cols, name = name) }

setMethod("show", "Bool", function(object) {
  size <- size(object)
  cat("Bool(", size[1], ", ", size[2], ")", sep = "")
})

setMethod("as.character", "Bool", function(x) {
  size <- size(x)
  paste("Bool(", size[1], ", ", size[2], ")", sep = "")
})

setMethod("canonicalize", "Bool", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constr <- canon[[2]]
  list(obj, c(constr, list(BoolConstr(obj))))
})

setMethod("is_positive", "Bool", function(object) { TRUE })
setMethod("is_negative", "Bool", function(object) { FALSE })

# Integer variable
.Int <- setClass("Int", contains = "Variable")
Int <- function(rows = 1, cols = 1, name = NA_character_) { .Int(rows = rows, cols = cols, name = name) }

setMethod("show", "Int", function(object) {
  size <- size(object)
  cat("Int(", size[1], ", ", size[2], ")", sep = "")
})

setMethod("as.character", "Int", function(x) {
  size <- size(x)
  paste("Int(", size[1], ", ", size[2], ")", sep = "")
})

setMethod("canonicalize", "Int", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constr <- canon[[2]]
  list(obj, c(constr, list(IntConstr(obj))))
})

# Non-negative variable
.NonNegative <- setClass("NonNegative", contains = "Variable")
NonNegative <- function(rows = 1, cols = 1, name = NA_character_) { .NonNegative(rows = rows, cols = cols, name = name) }

setMethod("show", "NonNegative", function(object) {
  size <- size(object)
  cat("NonNegative(", size[1], ", ", size[2], ")", sep = "")
})

setMethod("as.character", "NonNegative", function(x) {
  size <- size(x)
  paste("NonNegative(", size[1], ", ", size[2], ")", sep = "")
})

setMethod("canonicalize", "NonNegative", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constr <- canon[[2]]
  list(obj, c(constr, list(create_geq(obj))))
})

setMethod("is_positive", "NonNegative", function(object) { TRUE })
setMethod("is_negative", "NonNegative", function(object) { FALSE })

# Positive semidefinite matrix
.SemidefUpperTri <- setClass("SemidefUpperTri", representation(n = "numeric"), contains = "Variable")
SemidefUpperTri <- function(n, name = NA_character_) { .SemidefUpperTri(n = n, name = name) }

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

upper_tri_to_full <- function(n) {
  entries <- n*(n+1)/2
  
  val_arr <- c()
  row_arr <- c()
  col_arr <- c()
  count <- 1
  for(i in 1:n) {
    for(j in i:n) {
      # Index in the original matrix
      col_arr <- c(col_arr, count)
      
      # Index in the filled matrix
      row_arr <- c(row_arr, j*n + i)
      val_arr <- c(val_arr, 1.0)
      if(i != j) {
        # Index in the original matrix
        col_arr <- c(col_arr, count)
        
        # Index in the filled matrix
        row_arr <- c(row_arr, i*n + j)
        val_arr <- c(val_arr, 1.0)
      }
      count <- count + 1
    }
  }
  sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = c(n*n, entries))
}

setMethod("canonicalize", "SemidefUpperTri", function(object) {
  upper_tri <- create_var(c(size(object)[1], 1), object@id)
  fill_coef <- upper_tri_to_full(object@n)
  fill_coef <- create_const(fill_coef, c(object@n*object@n, size(object)[1]), sparse = TRUE)
  full_mat = mul_expr(fill_coef, upper_tri, c(object@n*object@n, 1))
  full_mat <- reshape(full_mat, c(object@n, object@n))
  list(upper_tri, list(SDP(full_mat, enforce_sym = FALSE)))
})

Semidef <- function(n, name) {
  var <- SemidefUpperTri(n, name)
  fill_mat <- Constant(upper_tri_to_full(n))
  Reshape(fill_mat*var, n, n)
}

# Symmetric matrix
.SymmetricUpperTri <- setClass("SymmetricUpperTri", representation(n = "numeric"), contains = "Variable")
SymmetricUpperTri <- function(n, name = NA_character_) { .SymmetricUpperTri(n = n, name = name) }

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

setMethod("canonicalize", "SymmetricUpperTri", function(object) {
  upper_tri <- create_var(c(size(object)[1], 1), object@id)
  list(upper_tri, list())
})

Symmetric <- function(n, name = NA_character_) {
  var <- SymmetricUpperTri(n, name)
  fill_mat <- Constant(upper_tri_to_full(n))
  Reshape(fill_mat*var, as.integer(n), as.integer(n))
}
