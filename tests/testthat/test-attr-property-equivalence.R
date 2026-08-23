## The invariant the `@`-bypass helpers rest on (ADR D_PERF.8).
##
## `.args()`, `.shape()`, `.id()` and `.attributes()` read a property with
## `attr(x, "<prop>", exact = TRUE)` instead of `x@<prop>`, because the S3
## dispatch behind `@` costs ~450ns per read and these four properties are 81.5%
## of all property reads on a solve.
##
## That is only correct while CVXR keeps three properties of its own:
##
##   1. every constructor goes through `.fast_new`, which writes properties as
##      PLAIN ATTRIBUTES (CLAUDE.md constraint 17);
##   2. no class has a computed property (a `getter=`);
##   3. no class has a property setter or validator that `@` would run.
##
## If any of those stops holding, the helpers silently return the raw stored
## value and skip the getter -- a wrong answer with NO error, in the hottest
## code in the package. Nothing in S7 or in the compiler catches it.
##
## So this file does not test four functions; it ENUMERATES every exported S7
## class in CVXR, constructs nothing, and checks the class metadata directly.
## It fails the day someone adds a getter, which is exactly when the helpers
## must be reverted.

## @cvxpy NONE
test_that("no CVXR S7 class has a computed property, setter or validator", {
  ns <- asNamespace("CVXR")
  classes <- Filter(function(nm) {
    obj <- get(nm, envir = ns)
    inherits(obj, "S7_class")
  }, ls(ns, all.names = TRUE))

  expect_gt(length(classes), 100L)   # the audit covered 162; guard the sweep

  offenders <- list()
  for (nm in classes) {
    cls <- get(nm, envir = ns)
    props <- S7::prop(cls, "properties")
    if (is.null(props)) next
    for (pn in names(props)) {
      p <- props[[pn]]
      bad <- character(0)
      if (!is.null(p$getter))    bad <- c(bad, "getter")
      if (!is.null(p$setter))    bad <- c(bad, "setter")
      if (!is.null(p$validator)) bad <- c(bad, "validator")
      if (length(bad))
        offenders[[length(offenders) + 1L]] <-
          sprintf("%s@%s has %s", nm, pn, paste(bad, collapse = "+"))
    }
    if (!is.null(S7::prop(cls, "validator")))
      offenders[[length(offenders) + 1L]] <- sprintf("%s has a class validator", nm)
  }

  expect_equal(unlist(offenders), NULL,
               info = paste0(
                 "A class gained a getter/setter/validator. The D_PERF.8 helpers ",
                 "(.args/.shape/.id/.attributes in zzz_R_specific/utility.R) read ",
                 "past it and MUST be reverted to `@` before this ships."))
})

## @cvxpy NONE
test_that("the helpers agree with `@` on live objects of many shapes", {
  x <- Variable(c(3L, 4L))
  p <- Parameter(2L, nonneg = TRUE)
  cst <- Constant(matrix(1:6, 2, 3))
  objs <- list(
    variable   = x,
    parameter  = p,
    constant   = cst,
    add        = x + x,
    negate     = -x,
    multiply   = 2 * x,
    sum        = sum_entries(x),
    index      = x[1L, 2L],
    transpose  = t(x),
    absolute   = abs(x),
    norm2      = p_norm(cst, 2),
    constraint = (x >= 0)
  )
  for (nm in names(objs)) {
    o <- objs[[nm]]
    if (!is.null(attr(o, "args", exact = TRUE)) || .s7_is(o, Expression))
      expect_identical(.args(o), o@args, info = nm)
    if (.s7_is(o, Expression)) {
      expect_identical(.shape(o), o@shape, info = nm)
      expect_identical(.id(o), o@id, info = nm)
    }
    if (.s7_is(o, Leaf))
      expect_identical(.attributes(o), o@attributes, info = nm)
  }
})

## @cvxpy NONE
test_that("exact = TRUE -- the helpers do not partial-match", {
  ## Without `exact = TRUE`, `attr(x, "arg")` returns `args`. A helper that
  ## partial-matched would be a silent wrong read, so pin the behavior that
  ## makes the helpers safe rather than trusting the argument is still there.
  x <- Variable(3L) + Variable(3L)
  expect_identical(attr(x, "arg", exact = TRUE), NULL)
  expect_identical(.args(x), x@args)
  expect_false(is.null(.args(x)))
})
