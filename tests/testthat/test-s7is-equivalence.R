## ## PERF-S7IS safety net -----------------------------------------------------
## Asserts the fast helpers .s7_is / .s7_is_any (rsrc_tree/zzz_R_specific/
## s7_dispatch_perf.R) are TAUTOLOGICALLY equal to S7::S7_inherits for every
## concrete/abstract CVXR class that the codebase tests membership against.
##
## This is the revert trigger: if S7 ever changes its class() layout (so that
## inherits(x, "CVXR::Cls") stops mirroring the S7 ancestor chain) this test
## FAILS LOUDLY in R CMD check, signaling that every `.s7_is` swap must be
## reverted to S7_inherits (one-line body change in s7_dispatch_perf.R).
## The class list below is a SUPERSET of all swapped sites (every distinct
## S7_inherits(., CLASS) literal in rsrc_tree), so it needs no upkeep as more
## sites are swapped.
## ----------------------------------------------------------------------------

ns <- asNamespace("CVXR")
.s7_is     <- get(".s7_is", ns)
.s7_is_any <- get(".s7_is_any", ns)

## --- panel of instances spanning the hierarchy --------------------------------
x       <- Variable(3)
p       <- Parameter()
cst     <- Constant(matrix(1:4, 2, 2))
add     <- x + 1                       # AddExpression
negx    <- -x                          # NegExpression
A       <- Constant(matrix(1, 2, 3))
mul     <- A %*% x                     # MulExpression
absx    <- abs(x)                      # Abs (Atom / non-AffAtom)
objmin  <- Minimize(sum(x))
objmax  <- Maximize(sum(-x))
prob    <- Problem(objmin, list(x >= 0))
con_eq  <- (x == 0)                     # Equality
con_ineq <- (x >= 0)                    # Inequality

instances <- list(
  x = x, p = p, cst = cst, add = add, negx = negx, mul = mul, absx = absx,
  objmin = objmin, objmax = objmax, prob = prob,
  con_eq = con_eq, con_ineq = con_ineq,
  num = 5, chr = "a", lst = list(1), nul = NULL
)

## AUTO-ENUMERATE every CVXR S7 class in the namespace. .s7_is swaps are no longer
## confined to a hand-maintained list -- after the full sweep, active code uses
## .s7_is everywhere -- so the safety net must cover EVERY CVXR class, present and
## future. This list therefore derives itself from the package, so it can never
## drift behind the swapped sites (a new class is covered the moment it is added).
classes <- Filter(Negate(is.null), lapply(ls(ns, all.names = TRUE), function(nm) {
  obj <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
  if (!is.null(obj) && inherits(obj, "S7_class") &&
      identical(attr(obj, "package"), "CVXR")) obj else NULL
}))
names(classes) <- vapply(classes, function(cls) attr(cls, "name"), character(1))

test_that(".s7_is == S7_inherits for every (instance, class) pair", {
  for (cls in classes) {
    for (inm in names(instances)) {
      inst <- instances[[inm]]
      expect_identical(
        .s7_is(inst, cls), S7::S7_inherits(inst, cls),
        info = sprintf("instance=%s class=%s", inm, attr(cls, "name"))
      )
    }
  }
})

test_that(".s7_is reproduces known positive memberships", {
  expect_true(.s7_is(x, classes[["Variable"]]))
  expect_true(.s7_is(x, classes[["Leaf"]]))
  expect_true(.s7_is(x, classes[["Expression"]]))
  expect_true(.s7_is(x, classes[["Canonical"]]))
  expect_true(.s7_is(cst, classes[["Constant"]]))
  expect_true(.s7_is(p, classes[["Parameter"]]))
  expect_true(.s7_is(add, classes[["AddExpression"]]))
  expect_true(.s7_is(prob, classes[["Problem"]]))
  expect_true(.s7_is(objmin, classes[["Minimize"]]))
  expect_true(.s7_is(objmax, classes[["Maximize"]]))
})

test_that(".s7_is is FALSE for non-S7 and unrelated inputs", {
  expect_false(.s7_is(NULL, classes[["Expression"]]))
  expect_false(.s7_is(5, classes[["Expression"]]))
  expect_false(.s7_is("a", classes[["Expression"]]))
  expect_false(.s7_is(list(1), classes[["Expression"]]))
  expect_false(.s7_is(x, classes[["Parameter"]]))   # Variable is not a Parameter
  expect_false(.s7_is(cst, classes[["Variable"]]))
})

test_that(".s7_is_any == any(S7_inherits) over class lists", {
  panel <- list(classes[["Variable"]], classes[["Parameter"]], classes[["Constant"]])
  for (inm in names(instances)) {
    inst <- instances[[inm]]
    expect_identical(
      .s7_is_any(inst, panel),
      any(vapply(panel, function(cls) S7::S7_inherits(inst, cls), logical(1))),
      info = sprintf("instance=%s", inm)
    )
  }
  expect_true(.s7_is_any(x, panel))
  expect_true(.s7_is_any(p, panel))
  expect_false(.s7_is_any(absx, panel))
  expect_false(.s7_is_any(NULL, panel))
})
