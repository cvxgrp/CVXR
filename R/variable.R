.Variable <- setClass("Variable", representation(rows = "numeric", cols = "numeric", name = "character"),
                                 prototype(rows = 1, cols = 1, name = NA_character_), contains = "Leaf")

Variable <- function(rows = 1, cols = 1, name = NA_character_) { .Variable(rows = rows, cols = cols, name = name) }

setMethod("init_dcp_attr", "Variable", function(object) {
  DCPAttr(sign = Sign(sign = SIGN_UNKNOWN_KEY), curvature = Curvature(curvature = CURV_AFFINE_KEY), shape = Shape(rows = object@rows, cols = object@cols))
})

setMethod("initialize", "Variable", function(.Object, ..., rows = 1, cols = 1, name = NA_character_) {
  .Object@rows = rows
  .Object@cols = cols
  .Object@name = name
  .Object@dcp_attr = init_dcp_attr(.Object)
  callNextMethod(.Object, ...)
})

setMethod("variables", "Variable", function(object) { list(object) })

setMethod("canonicalize", "Variable", function(object) {
  obj <- create_var(size(object))
  list(obj, list())
})

# Boolean variable
Bool <- setClass("Bool", contains = "Variable")

# Integer variable
Int <- setClass("Int", contains = "Variable")

# Non-negative variable
NonNegative <- setClass("NonNegative", contains = "Variable")

setMethod("init_dcp_attr", "NonNegative", function(object) {
  DCPAttr(Sign(sign = SIGN_POSITIVE_KEY), Curvature(curvature = CURV_AFFINE_KEY), Shape(rows = object@rows, cols = object@cols))
})

# Positive semidefinite matrix
SemidefUpperTri <- setClass("SemidefUpperTri", representation(n = "numeric"), contains = "Variable")

setMethod("initialize", "SemidefUpperTri", function(.Object, ..., rows, cols, name, n) {
  .Object@n = n
  callNextMethod(.Object, ..., rows = n*(n+1)/2, cols = 1, name = name)
})

# Symmetric matrix
SymmetricUpperTri <- setClass("SymmetricUpperTri", representation(n = "numeric"), contains = "Variable")

setMethod("initialize", "SemidefUpperTri", function(.Object, ..., rows, cols, name, n) {
  .Object@n = n
  callNextMethod(.Object, ..., rows = n*(n+1)/2, cols = 1, name = name)
})

