## Tests for Phase 2B: Atom + AffAtom base classes
##
## These tests verify:
## - Atom/AffAtom construction and validation
## - DCP composition rules (is_convex, is_concave)
## - Sign propagation (is_nonneg, is_nonpos)
## - Value evaluation via numeric_value
## - canonicalize for constant atoms
## - Complex argument handling
## - PSD/NSD propagation through AffAtom

# ── Test helpers: concrete Atom and AffAtom subclasses ──────────────
#
# Atom and AffAtom are abstract. We create minimal concrete subclasses
# for testing, implementing the required abstract generics.

## Helper to build an AffAtom-like object with new_object()
.make_affatom <- function(args, shape) {
  args <- lapply(args, as_expr)
  shape <- validate_shape(shape)
  obj <- new_object(S7_object(),
    id    = next_expr_id(),
    .cache = new.env(parent = emptyenv()),
    args  = args,
    shape = shape
  )
  validate_arguments(obj)
  obj
}

## Helper to build an Atom-like object with new_object()
.make_atom <- function(args, shape) {
  args <- lapply(args, as_expr)
  shape <- validate_shape(shape)
  obj <- new_object(S7_object(),
    id    = next_expr_id(),
    .cache = new.env(parent = emptyenv()),
    args  = args,
    shape = shape
  )
  validate_arguments(obj)
  obj
}

## A minimal AffAtom: identity function f(x) = x, shape = arg shape
## Increasing in arg 0, sign from args, affine.
TestIdentityAtom <- new_class("TestIdentityAtom", parent = AffAtom, package = "CVXR",
  constructor = function(x) {
    x <- as_expr(x)
    new_object(S7_object(),
      id    = next_expr_id(),
      .cache = new.env(parent = emptyenv()),
      args  = list(x),
      shape = validate_shape(x@shape)
    )
  }
)

method(shape_from_args, TestIdentityAtom) <- function(x) x@args[[1L]]@shape
method(numeric_value, TestIdentityAtom) <- function(x, values, ...) values[[1L]]
method(graph_implementation, TestIdentityAtom) <- function(x, arg_objs, shape, data = NULL, ...) {
  list(arg_objs[[1L]], list())
}

## A minimal AffAtom: negation f(x) = -x, shape = arg shape
## Decreasing in arg 0, flips sign.
TestNegAtom <- new_class("TestNegAtom", parent = AffAtom, package = "CVXR",
  constructor = function(x) {
    x <- as_expr(x)
    new_object(S7_object(),
      id    = next_expr_id(),
      .cache = new.env(parent = emptyenv()),
      args  = list(x),
      shape = validate_shape(x@shape)
    )
  }
)

method(shape_from_args, TestNegAtom) <- function(x) x@args[[1L]]@shape
method(is_incr, TestNegAtom) <- function(x, idx, ...) FALSE
method(is_decr, TestNegAtom) <- function(x, idx, ...) TRUE
method(sign_from_args, TestNegAtom) <- function(x) {
  c(is_nonneg = is_nonpos(x@args[[1L]]), is_nonpos = is_nonneg(x@args[[1L]]))
}
method(numeric_value, TestNegAtom) <- function(x, values, ...) -values[[1L]]
method(graph_implementation, TestNegAtom) <- function(x, arg_objs, shape, data = NULL, ...) {
  list(neg_expr_linop(arg_objs[[1L]]), list())
}

## A minimal convex (non-affine) atom: f(x) = x^2 (scalar only)
## Convex but NOT concave. Increasing for x >= 0, decreasing for x <= 0.
TestSquareAtom <- new_class("TestSquareAtom", parent = Atom, package = "CVXR",
  constructor = function(x) {
    x <- as_expr(x)
    obj <- new_object(S7_object(),
      id    = next_expr_id(),
      .cache = new.env(parent = emptyenv()),
      args  = list(x),
      shape = c(1L, 1L)
    )
    validate_arguments(obj)
    obj
  }
)

method(shape_from_args, TestSquareAtom) <- function(x) c(1L, 1L)
method(sign_from_args, TestSquareAtom) <- function(x) c(is_nonneg = TRUE, is_nonpos = FALSE)
method(is_atom_convex, TestSquareAtom) <- function(x) TRUE
method(is_atom_concave, TestSquareAtom) <- function(x) FALSE
method(is_incr, TestSquareAtom) <- function(x, idx, ...) is_nonneg(x@args[[1L]])
method(is_decr, TestSquareAtom) <- function(x, idx, ...) is_nonpos(x@args[[1L]])
method(numeric_value, TestSquareAtom) <- function(x, values, ...) values[[1L]]^2
method(graph_implementation, TestSquareAtom) <- function(x, arg_objs, shape, data = NULL, ...) {
  stop("graph_implementation not implemented for TestSquareAtom")
}

## A minimal concave (non-affine) atom: f(x) = -x^2 (scalar only)
TestNegSquareAtom <- new_class("TestNegSquareAtom", parent = Atom, package = "CVXR",
  constructor = function(x) {
    x <- as_expr(x)
    obj <- new_object(S7_object(),
      id    = next_expr_id(),
      .cache = new.env(parent = emptyenv()),
      args  = list(x),
      shape = c(1L, 1L)
    )
    validate_arguments(obj)
    obj
  }
)

method(shape_from_args, TestNegSquareAtom) <- function(x) c(1L, 1L)
method(sign_from_args, TestNegSquareAtom) <- function(x) c(is_nonneg = FALSE, is_nonpos = TRUE)
method(is_atom_convex, TestNegSquareAtom) <- function(x) FALSE
method(is_atom_concave, TestNegSquareAtom) <- function(x) TRUE
method(is_incr, TestNegSquareAtom) <- function(x, idx, ...) is_nonpos(x@args[[1L]])
method(is_decr, TestNegSquareAtom) <- function(x, idx, ...) is_nonneg(x@args[[1L]])
method(numeric_value, TestNegSquareAtom) <- function(x, values, ...) -(values[[1L]]^2)
method(graph_implementation, TestNegSquareAtom) <- function(x, arg_objs, shape, data = NULL, ...) {
  stop("graph_implementation not implemented for TestNegSquareAtom")
}

## A convex atom with NO monotonicity (neither incr nor decr).
## This makes DCP fail for non-affine args.
TestConvexNonMonotoneAtom <- new_class("TestConvexNonMonotoneAtom", parent = Atom, package = "CVXR",
  constructor = function(x) {
    x <- as_expr(x)
    obj <- new_object(S7_object(),
      id    = next_expr_id(),
      .cache = new.env(parent = emptyenv()),
      args  = list(x),
      shape = c(1L, 1L)
    )
    validate_arguments(obj)
    obj
  }
)

method(shape_from_args, TestConvexNonMonotoneAtom) <- function(x) c(1L, 1L)
method(sign_from_args, TestConvexNonMonotoneAtom) <- function(x) c(is_nonneg = FALSE, is_nonpos = FALSE)
method(is_atom_convex, TestConvexNonMonotoneAtom) <- function(x) TRUE
method(is_atom_concave, TestConvexNonMonotoneAtom) <- function(x) FALSE
method(is_incr, TestConvexNonMonotoneAtom) <- function(x, idx, ...) FALSE
method(is_decr, TestConvexNonMonotoneAtom) <- function(x, idx, ...) FALSE
method(numeric_value, TestConvexNonMonotoneAtom) <- function(x, values, ...) values[[1L]]
method(graph_implementation, TestConvexNonMonotoneAtom) <- function(x, arg_objs, shape, data = NULL, ...) {
  stop("graph_implementation not implemented for TestConvexNonMonotoneAtom")
}

## A two-argument AffAtom: f(x, y) = x + y (same shape)
TestAddAtom <- new_class("TestAddAtom", parent = AffAtom, package = "CVXR",
  constructor = function(x, y) {
    x <- as_expr(x)
    y <- as_expr(y)
    shape <- sum_shapes(list(x@shape, y@shape))
    new_object(S7_object(),
      id    = next_expr_id(),
      .cache = new.env(parent = emptyenv()),
      args  = list(x, y),
      shape = shape
    )
  }
)

method(shape_from_args, TestAddAtom) <- function(x) {
  sum_shapes(list(x@args[[1L]]@shape, x@args[[2L]]@shape))
}
method(numeric_value, TestAddAtom) <- function(x, values, ...) values[[1L]] + values[[2L]]
method(graph_implementation, TestAddAtom) <- function(x, arg_objs, shape, data = NULL, ...) {
  list(sum_expr_linop(arg_objs), list())
}


# ═══════════════════════════════════════════════════════════════════
# Tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Atom class exists and has correct inheritance", {
  expect_true(inherits(Atom, "S7_class"))
  ## Atom → Expression → Canonical
  x <- Variable(1)
  a <- TestIdentityAtom(x)
  expect_true(S7_inherits(a, Atom))
  expect_true(S7_inherits(a, Expression))
  expect_true(S7_inherits(a, Canonical))
})

## @cvxpy NONE
test_that("AffAtom class exists and has correct inheritance", {
  expect_true(inherits(AffAtom, "S7_class"))
  x <- Variable(1)
  a <- TestIdentityAtom(x)
  expect_true(S7_inherits(a, AffAtom))
  expect_true(S7_inherits(a, Atom))
})

## @cvxpy NONE
test_that("Atom constructor rejects zero arguments", {
  expect_error(Atom(args = list(), shape = c(1L, 1L)),
               "No arguments given")
})

## @cvxpy NONE
test_that("AffAtom constructor rejects zero arguments", {
  expect_error(AffAtom(args = list(), shape = c(1L, 1L)),
               "No arguments given")
})

## @cvxpy NONE
test_that("Atom constructor converts raw values to Constant", {
  ## Pass a plain numeric — should be converted via as_expr
  a <- TestIdentityAtom(5.0)
  expect_true(S7_inherits(a@args[[1L]], Constant))
  expect_equal(drop(value(a@args[[1L]])), 5.0)
})

## @cvxpy NONE
test_that("Atom constructor accepts Expression arguments", {
  x <- Variable(3)
  a <- TestIdentityAtom(x)
  expect_identical(a@args[[1L]], x)
})

## @cvxpy NONE
test_that("Atom has unique ID", {
  a1 <- TestIdentityAtom(1)
  a2 <- TestIdentityAtom(2)
  expect_true(a1@id != a2@id)
})

## @cvxpy NONE
test_that("Atom has correct shape", {
  x <- Variable(c(3, 4))
  a <- TestIdentityAtom(x)
  expect_equal(a@shape, c(3L, 4L))
})

## @cvxpy NONE
test_that("Atom has cache environment", {
  a <- TestIdentityAtom(1)
  expect_true(is.environment(a@.cache))
})

# ── validate_arguments tests ────────────────────────────────────────

## @cvxpy NONE
test_that("Atom.validate_arguments rejects complex args", {
  x <- Variable(1, complex = TRUE)
  expect_error(TestSquareAtom(x), "cannot be complex")
})

## @cvxpy NONE
test_that("AffAtom.validate_arguments allows complex args", {
  x <- Variable(1, complex = TRUE)
  a <- TestIdentityAtom(x)
  expect_true(S7_inherits(a, AffAtom))
})

# ── expr_name tests ─────────────────────────────────────────────────

## @cvxpy NONE
test_that("Atom expr_name includes class name and arg names", {
  x <- Variable(1, name = "x")
  a <- TestIdentityAtom(x)
  name <- expr_name(a)
  expect_true(grepl("TestIdentityAtom", name))
  expect_true(grepl("x", name))
})

# ── Sign propagation tests ──────────────────────────────────────────

## @cvxpy NONE
test_that("AffAtom sign propagation: positive args", {
  x <- Variable(1, nonneg = TRUE)
  a <- TestIdentityAtom(x)
  expect_true(is_nonneg(a))
  expect_false(is_nonpos(a))
})

## @cvxpy NONE
test_that("AffAtom sign propagation: negative args", {
  x <- Variable(1, nonpos = TRUE)
  a <- TestIdentityAtom(x)
  expect_false(is_nonneg(a))
  expect_true(is_nonpos(a))
})

## @cvxpy NONE
test_that("AffAtom sign propagation: zero args", {
  ## Constant 0 is both nonneg and nonpos
  a <- TestIdentityAtom(0)
  expect_true(is_nonneg(a))
  expect_true(is_nonpos(a))
  expect_true(is_zero(a))
})

## @cvxpy NONE
test_that("AffAtom sign propagation: unknown sign", {
  x <- Variable(1)
  a <- TestIdentityAtom(x)
  expect_false(is_nonneg(a))
  expect_false(is_nonpos(a))
})

## @cvxpy NONE
test_that("TestNegAtom flips sign", {
  x <- Variable(1, nonneg = TRUE)
  a <- TestNegAtom(x)
  expect_false(is_nonneg(a))
  expect_true(is_nonpos(a))
})

## @cvxpy NONE
test_that("Sign from constant atom", {
  a <- TestSquareAtom(Constant(3.0))
  expect_true(is_nonneg(a))
  expect_false(is_nonpos(a))
})

# ── Curvature: AffAtom ──────────────────────────────────────────────

## @cvxpy NONE
test_that("AffAtom is_atom_convex and is_atom_concave", {
  x <- Variable(1)
  a <- TestIdentityAtom(x)
  expect_true(is_atom_convex(a))
  expect_true(is_atom_concave(a))
  expect_true(is_atom_affine(a))
})

## @cvxpy NONE
test_that("AffAtom of variable is affine", {
  x <- Variable(1)
  a <- TestIdentityAtom(x)
  expect_true(is_convex(a))
  expect_true(is_concave(a))
  expect_true(is_affine(a))
  expect_true(is_dcp(a))
})

## @cvxpy NONE
test_that("AffAtom of constant is constant", {
  a <- TestIdentityAtom(5.0)
  expect_true(is_constant(a))
  expect_true(is_convex(a))
  expect_true(is_concave(a))
  expect_true(is_affine(a))
})

# ── DCP composition: convex atom ────────────────────────────────────

## @cvxpy NONE
test_that("Convex atom of affine arg is convex", {
  x <- Variable(1)
  a <- TestSquareAtom(x)
  expect_true(is_convex(a))
  expect_false(is_concave(a))
  expect_true(is_dcp(a))
})

## @cvxpy NONE
test_that("Convex atom of convex arg with incr monotonicity is convex", {
  ## square(x) where x is nonneg variable → increasing
  x <- Variable(1, nonneg = TRUE)
  inner <- TestSquareAtom(x)  ## convex, nonneg
  ## square(inner) = square(square(x)) = x^4
  ## is_incr checks if inner is nonneg → yes
  outer <- TestSquareAtom(inner)
  expect_true(is_convex(outer))
})

## @cvxpy NONE
test_that("Convex atom of concave arg with decr monotonicity IS convex", {
  ## square(neg_square(x)) — concave arg, square IS decr (arg is nonpos) → convex
  x <- Variable(1)
  inner <- TestNegSquareAtom(x)  ## concave, nonpositive
  outer <- TestSquareAtom(inner)
  ## square is decreasing for nonpos args → convex(concave + decr) → convex
  expect_true(is_convex(outer))
  expect_true(is_dcp(outer))
})

## @cvxpy NONE
test_that("Non-monotone convex atom of convex arg is NOT DCP", {
  ## ConvexNonMonotone(square(x)) — convex arg, not incr/decr → fails DCP
  x <- Variable(1)
  inner <- TestSquareAtom(x)   ## convex, not affine
  outer <- TestConvexNonMonotoneAtom(inner)
  expect_false(is_convex(outer))
  expect_false(is_dcp(outer))
})

## @cvxpy NONE
test_that("Convex atom of constant is convex", {
  a <- TestSquareAtom(Constant(3.0))
  expect_true(is_convex(a))
  expect_true(is_concave(a))  ## constant → both
  expect_true(is_constant(a))
})

# ── DCP composition: concave atom ──────────────────────────────────

## @cvxpy NONE
test_that("Concave atom of affine arg is concave", {
  x <- Variable(1)
  a <- TestNegSquareAtom(x)
  expect_false(is_convex(a))
  expect_true(is_concave(a))
  expect_true(is_dcp(a))
})

## @cvxpy NONE
test_that("Concave atom of concave arg with incr monotonicity is concave", {
  ## neg_square(x) where x is nonpos variable → is_incr
  x <- Variable(1, nonpos = TRUE)
  inner <- TestNegSquareAtom(x)  ## concave, nonpos
  ## neg_square(inner) — is_incr checks if inner is nonpos → yes
  outer <- TestNegSquareAtom(inner)
  expect_true(is_concave(outer))
})

# ── DCP composition: affine atom chain ──────────────────────────────

## @cvxpy NONE
test_that("Affine chain: identity(identity(x)) is affine", {
  x <- Variable(3)
  a <- TestIdentityAtom(TestIdentityAtom(x))
  expect_true(is_affine(a))
  expect_true(is_dcp(a))
})

## @cvxpy NONE
test_that("Affine chain: neg(neg(x)) is affine", {
  x <- Variable(3)
  a <- TestNegAtom(TestNegAtom(x))
  expect_true(is_affine(a))
})

# ── Two-argument atom tests ─────────────────────────────────────────

## @cvxpy NONE
test_that("TestAddAtom shape is correct", {
  x <- Variable(3)
  y <- Variable(3)
  a <- TestAddAtom(x, y)
  expect_equal(a@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("TestAddAtom is affine", {
  x <- Variable(3)
  y <- Variable(3)
  a <- TestAddAtom(x, y)
  expect_true(is_affine(a))
  expect_true(is_dcp(a))
})

## @cvxpy NONE
test_that("TestAddAtom sign: both nonneg", {
  x <- Variable(1, nonneg = TRUE)
  y <- Variable(1, nonneg = TRUE)
  a <- TestAddAtom(x, y)
  expect_true(is_nonneg(a))
  expect_false(is_nonpos(a))
})

## @cvxpy NONE
test_that("TestAddAtom sign: mixed signs → unknown", {
  x <- Variable(1, nonneg = TRUE)
  y <- Variable(1, nonpos = TRUE)
  a <- TestAddAtom(x, y)
  expect_false(is_nonneg(a))
  expect_false(is_nonpos(a))
})

# ── value evaluation ────────────────────────────────────────────────

## @cvxpy NONE
test_that("Atom value: constant args", {
  a <- TestIdentityAtom(5.0)
  expect_equal(drop(value(a)), 5.0)
})

## @cvxpy NONE
test_that("Atom value: square of constant", {
  a <- TestSquareAtom(3.0)
  expect_equal(drop(value(a)), 9.0)
})

## @cvxpy NONE
test_that("Atom value: neg_square of constant", {
  a <- TestNegSquareAtom(3.0)
  expect_equal(drop(value(a)), -9.0)
})

## @cvxpy NONE
test_that("Atom value: add two constants", {
  a <- TestAddAtom(3.0, 4.0)
  expect_equal(drop(value(a)), 7.0)
})

## @cvxpy NONE
test_that("Atom value: variable without value → NULL", {
  x <- Variable(1)
  a <- TestIdentityAtom(x)
  expect_null(value(a))
})

## @cvxpy NONE
test_that("Atom value: variable with value → evaluates", {
  x <- Variable(1)
  value(x) <- 5.0
  a <- TestIdentityAtom(x)
  expect_equal(drop(value(a)), 5.0)
})

## @cvxpy NONE
test_that("Atom value: nested atoms evaluate correctly", {
  ## square(constant(3)) = 9
  a <- TestSquareAtom(Constant(3.0))
  expect_equal(drop(value(a)), 9.0)
  ## identity(square(constant(3))) = 9
  b <- TestIdentityAtom(a)
  expect_equal(drop(value(b)), 9.0)
})

## @cvxpy NONE
test_that("Atom value: parameter without value → NULL", {
  p <- Parameter(1)
  a <- TestIdentityAtom(p)
  expect_null(value(a))
})

## @cvxpy NONE
test_that("Atom value: parameter with value → evaluates", {
  p <- Parameter(1)
  value(p) <- 7.0
  a <- TestIdentityAtom(p)
  expect_equal(drop(value(a)), 7.0)
})

# ── canonicalize tests ──────────────────────────────────────────────

## @cvxpy NONE
test_that("Atom canonicalize: constant atom wraps as Constant", {
  a <- TestIdentityAtom(5.0)
  cf <- canonicalize(a)
  expect_true(is.list(cf))
  expect_equal(length(cf), 2L)
  ## First element is a LinOp, second is constraints (list)
  expect_true(is.list(cf[[2L]]))
})

## @cvxpy NONE
test_that("Atom canonical_form is cached", {
  a <- TestIdentityAtom(5.0)
  cf1 <- canonical_form(a)
  cf2 <- canonical_form(a)
  ## Same result (cached)
  expect_identical(cf1, cf2)
})

# ── atoms() collection tests ────────────────────────────────────────

## @cvxpy NONE
test_that("Atom.atoms() returns class of the atom", {
  x <- Variable(1)
  a <- TestIdentityAtom(x)
  at <- atoms(a)
  expect_true(length(at) >= 1L)
  ## Should contain TestIdentityAtom class
  names <- vapply(at, function(cls) {
    if (inherits(cls, "S7_class")) cls@name else as.character(cls)
  }, character(1))
  expect_true("TestIdentityAtom" %in% names)
})

## @cvxpy NONE
test_that("Atom.atoms() collects from nested atoms", {
  x <- Variable(1)
  inner <- TestIdentityAtom(x)
  outer <- TestNegAtom(inner)
  at <- atoms(outer)
  names <- vapply(at, function(cls) {
    if (inherits(cls, "S7_class")) cls@name else as.character(cls)
  }, character(1))
  expect_true("TestNegAtom" %in% names)
  expect_true("TestIdentityAtom" %in% names)
})

## @cvxpy NONE
test_that("Leaf.atoms() returns empty list", {
  x <- Variable(1)
  expect_equal(length(atoms(x)), 0L)
})

# ── domain tests ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Atom.domain() returns empty list for basic atoms", {
  x <- Variable(1)
  a <- TestIdentityAtom(x)
  expect_equal(domain(a), list())
})

# ── variables / parameters / constants collection ───────────────────

## @cvxpy NONE
test_that("Atom collects variables from args", {
  x <- Variable(3)
  y <- Variable(3)
  a <- TestAddAtom(x, y)
  vars <- variables(a)
  expect_equal(length(vars), 2L)
})

## @cvxpy NONE
test_that("Atom collects parameters from args", {
  p <- Parameter(1)
  x <- Variable(1)
  a <- TestAddAtom(p, x)
  params <- parameters(a)
  expect_equal(length(params), 1L)
})

## @cvxpy NONE
test_that("Atom collects constants from args", {
  a <- TestAddAtom(Constant(1), Constant(2))
  consts <- constants(a)
  expect_equal(length(consts), 2L)
})

# ── PSD/NSD propagation through AffAtom ─────────────────────────────

## @cvxpy NONE
test_that("AffAtom PSD propagates: identity of PSD is PSD", {
  x <- Variable(c(3, 3), PSD = TRUE)
  a <- TestIdentityAtom(x)
  expect_true(is_psd(a))
  expect_false(is_nsd(a))
})

## @cvxpy NONE
test_that("AffAtom NSD propagates: identity of NSD is NSD", {
  x <- Variable(c(3, 3), NSD = TRUE)
  a <- TestIdentityAtom(x)
  expect_false(is_psd(a))
  expect_true(is_nsd(a))
})

## @cvxpy NONE
test_that("AffAtom PSD through negation: neg of NSD is PSD", {
  x <- Variable(c(3, 3), NSD = TRUE)
  a <- TestNegAtom(x)
  ## TestNegAtom is decreasing → is_psd needs decr(idx) && arg.is_nsd()
  ## But is_psd on AffAtom checks: incr(idx)&&psd OR decr(idx)&&nsd
  ## TestNegAtom: is_decr → TRUE, arg is NSD → TRUE
  ## Wait: is_psd checks decr(idx)&&arg.is_nsd() → that should be for NSD
  ## Actually: is_psd: (incr and psd) or (decr and nsd)
  ## decr and nsd → nsd check... no.
  ## Let me re-read:
  ## AffAtom.is_psd: for each arg: (incr(idx) && arg.is_psd()) || (decr(idx) && arg.is_nsd())
  ## TestNegAtom: decr(0) = TRUE, arg is NSD
  ## So decr(0) && arg.is_nsd() → TRUE → is_psd = TRUE
  expect_true(is_psd(a))
})

## @cvxpy NONE
test_that("AffAtom NSD through negation: neg of PSD is NSD", {
  x <- Variable(c(3, 3), PSD = TRUE)
  a <- TestNegAtom(x)
  ## AffAtom.is_nsd: (decr(idx) && arg.is_psd()) || (incr(idx) && arg.is_nsd())
  ## TestNegAtom: decr(0) = TRUE, arg is PSD → TRUE
  expect_true(is_nsd(a))
})

## @cvxpy NONE
test_that("AffAtom PSD: non-PSD/NSD arg → not PSD", {
  x <- Variable(c(3, 3))
  a <- TestIdentityAtom(x)
  expect_false(is_psd(a))
  expect_false(is_nsd(a))
})

# ── Complex propagation through AffAtom ─────────────────────────────

## @cvxpy NONE
test_that("AffAtom is_complex: complex arg → TRUE", {
  x <- Variable(1, complex = TRUE)
  a <- TestIdentityAtom(x)
  expect_true(is_complex(a))
})

## @cvxpy NONE
test_that("AffAtom is_complex: real arg → FALSE", {
  x <- Variable(1)
  a <- TestIdentityAtom(x)
  expect_false(is_complex(a))
})

## @cvxpy NONE
test_that("AffAtom is_imag: imaginary arg → TRUE", {
  x <- Variable(1, imag = TRUE)
  a <- TestIdentityAtom(x)
  expect_true(is_imag(a))
})

## @cvxpy NONE
test_that("AffAtom is_imag: mixed real/imag → FALSE", {
  x <- Variable(1, imag = TRUE)
  y <- Variable(1)
  a <- TestAddAtom(x, y)
  expect_false(is_imag(a))
  expect_true(is_complex(a))  ## any complex → TRUE
})

# ── is_atom_affine helper function ──────────────────────────────────

## @cvxpy NONE
test_that("is_atom_affine works", {
  x <- Variable(1)
  aff <- TestIdentityAtom(x)
  expect_true(is_atom_affine(aff))
  cvx <- TestSquareAtom(x)
  expect_false(is_atom_affine(cvx))
})

# ── Edge cases ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Atom with matrix argument", {
  m <- matrix(1:6, nrow = 3, ncol = 2)
  a <- TestIdentityAtom(m)
  expect_equal(a@shape, c(3L, 2L))
  expect_equal(as.matrix(value(a)), m)
})

## @cvxpy NONE
test_that("Atom sign is cached", {
  x <- Variable(1, nonneg = TRUE)
  a <- TestIdentityAtom(x)
  r1 <- is_nonneg(a)
  r2 <- is_nonneg(a)
  expect_identical(r1, r2)
  expect_true(cache_has(a, "is_nonneg"))
})

## @cvxpy NONE
test_that("Atom curvature is cached", {
  x <- Variable(1)
  a <- TestSquareAtom(x)
  r1 <- is_convex(a)
  r2 <- is_convex(a)
  expect_identical(r1, r2)
  expect_true(cache_has(a, "is_convex"))
})

## @cvxpy NONE
test_that("Atom print works without error", {
  x <- Variable(1)
  a <- TestIdentityAtom(x)
  expect_output(print(a), "TestIdentityAtom")
})
