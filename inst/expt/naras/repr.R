setMethod("repr", "Atom", function(object) {
    sprintf("%s(%s)", class(object),
            paste(lapply(object@args, as.character), collapse = ", "))
})

setMethod("repr", "Constant", function(object) {
    sprintf("Constant(%s, %s, (%s))",
            curvature(object), sign(object), paste(size(object), collapse = ","))
})

setMethod("repr", "Parameter", function(object) {
    sprintf('Parameter(%d, %d, sign="%s")',
            object@rows, object@cols, sign(object))
})

setMethod("repr", "LeqConstraint", function(object) {
    sprintf("%s(%s, %s)", class(object), repr(object@args[[1]]), repr(object@args[2]))
})

setMethod("repr", "Expression", function(object) {
    sprintf("Expression(%s, %s, %s)", curvature(object), sign(object),
            paste(size(object), collapse = ","))
})

setMethod("repr", "Variable", function(object) {
    size <- size(object)
    sprintf("Variable(%d, %d)", size[1], size[2])
})

setMethod("repr", "Bool", function(object) {
    size <- size(object)
    sprintf("Bool(%d, %d)", size[1], size[2])
})

setMethod("repr", "Int", function(object) {
    size <- size(object)
    sprintf("Int(%d, %d)", size[1], size[2])
})

setMethod("repr", "NonNegative", function(object) {
    size <- size(object)
    sprintf("NonNegative(%d, %d)", size[1], size[2])
})

setMethod("repr", "SemidefUpperTri", function(object) {
    sprintf("SemidefUpperTri(%d)", object@n)
})

setMethod("repr", "SymmetricUpperTri", function(object) {
    sprintf("SymmetricUpperTri(%d)", object@n)
})


