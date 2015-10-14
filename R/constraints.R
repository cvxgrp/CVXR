setClass("Constraint", representation(id = "character"), prototype(id = NA_character_), contains = "VIRTUAL")
setMethod("id", "Constraint", function(object) { object@id })

# TODO: What type should lin_op be? Need linear operator object.
BoolConstr <- setClass("BoolConstr", representation(lin_op = "ANY", .noncvx_var = "LinOp"), contains = "Constraint")

setMethod("initialize", "BoolConstr", function(.Object, ..., lin_op, .noncvx_var) {
  .Object@lin_op <- lin_op
  if(lin_op == VARIABLE)
    .Object@.noncvx_var <- .Object@lin_op
  else
    .Object@.noncvx_var <- create_var(size(.Object@lin_op))
  callNextMethod(.Object, ...)
})

setMethod("size", "BoolConstr", function(object) { size(.Object@lin_op) })

IntConstr <- setClass("IntConstr", contains = "BoolConstr")

LeqConstraint <- setClass("LeqConstraint", representation(lh_exp = "Expression", rh_exp = "Expression", .expr = "Expression"),
                           prototype(lh_exp = new("Expression"), rh_exp = new("Expression"), .expr = NULL),
                           validity = function(object) {
                             if(!is.null(object@.expr))
                               stop("[LeqConstraint: .expr] .expr is an internal slot and should not be set by user")
                             return(TRUE)
                           },
                            contains = "Constraint")

setMethod("initialize", "LeqConstraint", definition = function(.Object, id, lh_exp, rh_exp, .expr) {
  .Object@id <- id
  .Object@lh_exp <- lh_exp
  .Object@rh_exp <- rh_exp
  .Object@.expr <- lh_exp - rh_exp
  return(.Object)
})

setMethod("size", "LeqConstraint", function(object) { size(object@.expr) })

setMethod("is_dcp", "LeqConstraint", function(object) { is_convex(object@.expr) })

setMethod("canonicalize", "LeqConstraint", function(object) {
  canon <- canonical_form(object@.expr)
  dual <- create_leq(canon[[1]], constr_id = object@id)
  c(NA, list(canon[[2]], dual))
})

setMethod("variables", "LeqConstraint", function(object) { variables(object@.expr) })

setMethod("parameters", "LeqConstraint", function(object) { parameters(object@.expr) })

EqConstraint <- setClass("EqConstraint", contains = "LeqConstraint")

setMethod("is_dcp", "EqConstraint", function(object) { is_affine(object@.expr) })

setMethod("canonicalize", "EqConstraint", function(object) {
  canon <- canonical_form(object@.expr)
  dual <- create_eq(canon[[1]], constr_id = object@id)
  c(NA, list(canon[[2]], dual))
})

PSDConstraint <- setClass("PSDConstraint", contains = "LeqConstraint", 
                           validity = function(object) {
                             lh_exp <- object@lh_exp
                             rh_exp <- object@rh_exp
                             if(lh_exp@size[1] != lh_exp@size[2] || rh_exp@size[1] != rh_exp@size[2])
                               stop("[PSDConstraint: validation] non-square matrix in positive definite constraint")
                             return(TRUE)
                           })

setMethod("is_dcp", "PSDConstraint", function(object) { is_affine(object@.expr) })

SOC <- setClass("SOC", representation(t = "numeric", x_elems = "numeric"), 
                       prototype(t = NA_real_, x_elems = NA_real_), contains = "Constraint")

setMethod("show", "SOC", function(object) { cat("SOC(", object@t, ", <", paste(object@x_elems, collapse = ","), ">)", sep = "") })

setMethod("size", "SOC", function(object) {
  sizes <- sapply(object@x_elems, function(elem) { size(elem)[1] * size(elem)[2] })
  rows <- sum(sizes) + 1
  c(rows, 1)
})

SDP <- setClass("SDP", representation(A = "numeric", enforce_sym = "logical"),
                       prototype(A = NA_real_, enforce_sym = TRUE), contains = "Constraint")

setMethod("show", "SDP", function(object) { cat("SDP(", object@A, ")", sep = "") })

setMethod("size", "SDP", function(object) { size(object@A) })

SOCElemwise <- setClass("SOCElemwise", contains = "SOC")

setMethod("num_cones", "SOCElemwise", function(object) {
  size(object@t)[1] * size(object@t)[2]
})

setMethod("cone_size", "SOCElemwise", function(object) {
  c(1 + length(object@x_elems), 1)
})

setMethod("size", "SOCElemwise", function(object) {
  cone_size <- cone_size(object)
  rep(cone_size, num_cones(object))
})
