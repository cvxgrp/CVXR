.Variable <- setClass("Variable", representation(id = "character", rows = "numeric", cols = "numeric", name = "character", primal_value = "numeric"),
                                 prototype(id = UUIDgenerate(), rows = 1, cols = 1, name = NA_character_, primal_value = NA_real_), 
                                 validity = function(object) {
                                   if(!is.na(object@primal_value))
                                     stop("[Variable: primal_value] primal_value is an internal slot and should not be set by user")
                                   return(TRUE)
                                 }, contains = "Leaf")

Variable <- function(rows = 1, cols = 1, name = NA_character_) { .Variable(rows = rows, cols = cols, name = name) }

setMethod("init_dcp_attr", "Variable", function(object) {
  DCPAttr(sign = Sign(sign = SIGN_UNKNOWN_KEY), curvature = Curvature(curvature = CURV_AFFINE_KEY), shape = Shape(rows = object@rows, cols = object@cols))
})

setMethod("initialize", "Variable", function(.Object, ..., rows = 1, cols = 1, name = NA_character_) {
  .Object@rows <- rows
  .Object@cols <- cols
  if(is.na(name))
    .Object@name <- sprintf("%s%s", VAR_PREFIX, .Object@id)
  else
    .Object@name <- name
  .Object@dcp_attr = init_dcp_attr(.Object)
  callNextMethod(.Object, ...)
})

setMethod("get_data", "Variable", function(object) { list(object@rows, object@cols, object@name) })
setMethod("name", "Variable", function(object) { object@name })
setMethod("variables", "Variable", function(object) { list(object) })
setMethod("value", "Variable", function(object) { object@primal_value })

setMethod("save_value", "Variable", function(object) { 
  object@primal_value <- value
  object
})

setMethod("canonicalize", "Variable", function(object) {
  obj <- create_var(size(object))
  list(obj, list())
})

# Boolean variable
Bool <- setClass("Bool", contains = "Variable")
setMethod("canonicalize", "Bool", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constr <- canon[[2]]
  list(obj, c(constr, BoolConstr(obj)))
})

# Integer variable
Int <- setClass("Int", contains = "Variable")
setMethod("canonicalize", "Int", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constr <- canon[[2]]
  list(obj, c(constr, IntConstr(obj)))
})

# Non-negative variable
NonNegative <- setClass("NonNegative", contains = "Variable")

setMethod("init_dcp_attr", "NonNegative", function(object) {
  DCPAttr(Sign(sign = SIGN_POSITIVE_KEY), Curvature(curvature = CURV_AFFINE_KEY), Shape(rows = object@rows, cols = object@cols))
})

setMethod("canonicalize", "NonNegative", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constr <- canon[[2]]
  list(obj, c(constr, create_geq(obj)))
})

# Positive semidefinite matrix
.SemidefUpperTri <- setClass("SemidefUpperTri", representation(n = "numeric"), contains = "Variable")
SemidefUpperTri <- function(n, name) { 
  if(missing(name))
    .SemidefUpperTri(n = n) 
  else
    .SemidefUpperTri(n = n, name = name)
}

setMethod("initialize", "SemidefUpperTri", function(.Object, ..., rows, cols, n) {
  .Object@n = n
  callNextMethod(.Object, ..., rows = n*(n+1)/2, cols = 1)
})

setMethod("canonicalize", "SemidefUpperTri", function(object) {
  upper_tri <- create_var(c(size(object)[1], 1), object@id)
  fill_coef <- upper_tri_to_full(object@n)
  fill_coef <- create_const(fill_coef, c(object@n*object@n, size(object)[1]))
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
SymmetricUpperTri <- setClass("SymmetricUpperTri", representation(n = "numeric"), contains = "Variable")

setMethod("initialize", "SemidefUpperTri", function(.Object, ..., rows, cols, n) {
  .Object@n = n
  callNextMethod(.Object, ..., rows = n*(n+1)/2, cols = 1)
})

setMethod("canonicalize", "SemidefUpperTri", function(object) {
  upper_tri <- create_var(c(size(object)[1], 1), object@id)
  list(upper_tri, list())
})

Symmetric <- function(n, name) {
  var <- SymmetricUpperTri(n, name)
  fill_mat <- Constant(upper_tri_to_full(n))
  Reshape(fill_mat*var, round(n), round(n))
}
