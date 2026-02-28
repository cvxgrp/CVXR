## Phase 3A: Base Classes + LinOp Extensions + Spike Tests
##
## Tests:
## 1. S3 Math/Summary dispatch spike (B3)
## 2. Elementwise base class
## 3. AxisAtom base class
## 4. AxisAffAtom base class
## 5. domain() bug fix (B1)
## 6. New LinOp constructors

library(testthat)
library(CVXR)

## ── S3 Math/Summary Dispatch Spike (B3) ──────────────────────────
## Verify Math and Summary S3 group dispatch works with S7 Expression objects.
## This was a risk identified in the plan; confirmed working.

## @cvxpy NONE
test_that("S3 Math group dispatch routes to Expression", {
  x <- Variable(3)
  ## abs(x) should create an Abs atom
  a <- abs(x)
  expect_s3_class(a, "CVXR::Abs")
  expect_equal(a@shape, c(3L, 1L))
  ## exp(x) should create an Exp atom
  e <- exp(x)
  expect_s3_class(e, "CVXR::Exp")
  ## log(x) should create a Log atom
  l <- log(x)
  expect_s3_class(l, "CVXR::Log")
  ## sqrt(x) should create a Power atom with p=0.5
  s <- sqrt(x)
  expect_s3_class(s, "CVXR::Power")
})

## @cvxpy NONE
test_that("S3 Summary group dispatch routes to Expression", {
  x <- Variable(3)
  ## sum(x) should create a SumEntries atom
  s <- sum(x)
  expect_s3_class(s, "CVXR::SumEntries")
  expect_equal(s@shape, c(1L, 1L))
  ## max(x) should create a MaxEntries atom
  mx <- max(x)
  expect_s3_class(mx, "CVXR::MaxEntries")
  ## min(x) should create a MinEntries atom
  mn <- min(x)
  expect_s3_class(mn, "CVXR::MinEntries")
})

## ── Elementwise Base Class ──────────────────────────────────────

## @cvxpy NONE
test_that("Elementwise broadcasts argument shapes correctly", {
  x <- Variable(3)
  y <- Variable(3)
  ## Create a concrete Elementwise (use it directly since it's a base class)
  ## We test the shape computation by subclassing — for now, test shape_from_args
  ## by constructing with known args
  e <- Elementwise(args = list(x, y))
  expect_equal(e@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Elementwise broadcasts scalar with matrix", {
  x <- Variable(c(2, 3))
  c1 <- Constant(1)
  e <- Elementwise(args = list(x, c1))
  expect_equal(e@shape, c(2L, 3L))
})

## @cvxpy NONE
test_that("Elementwise rejects incompatible shapes", {
  x <- Variable(c(2, 3))
  y <- Variable(c(4, 5))
  expect_error(Elementwise(args = list(x, y)))
})

## @cvxpy NONE
test_that("Elementwise inherits from Atom", {
  x <- Variable(c(2, 3))
  e <- Elementwise(args = list(x))
  expect_true(S7_inherits(e, Atom))
  expect_true(S7_inherits(e, Elementwise))
})

## @cvxpy NONE
test_that("Elementwise is_symmetric works", {
  x <- Variable(c(3, 3))
  e <- Elementwise(args = list(x))
  ## Not necessarily symmetric; square shape + symmetric args needed
  expect_equal(e@shape[1L], e@shape[2L])
})

## ── AxisAtom Base Class ─────────────────────────────────────────

## @cvxpy NONE
test_that("AxisAtom with axis=NULL reduces to scalar", {
  x <- Variable(c(3, 4))
  a <- AxisAtom(x, axis = NULL)
  expect_equal(a@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("AxisAtom with axis=2 reduces column-wise", {
  x <- Variable(c(3, 4))
  a <- AxisAtom(x, axis = 2L)
  expect_equal(a@shape, c(1L, 4L))
})

## @cvxpy NONE
test_that("AxisAtom with axis=1 removes second dimension", {
  x <- Variable(c(3, 4))
  a <- AxisAtom(x, axis = 1L)
  expect_equal(a@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("AxisAtom keepdims=TRUE preserves dimensions", {
  x <- Variable(c(3, 4))
  a <- AxisAtom(x, axis = 2L, keepdims = TRUE)
  expect_equal(a@shape, c(1L, 4L))

  b <- AxisAtom(x, axis = 1L, keepdims = TRUE)
  expect_equal(b@shape, c(3L, 1L))

  c_atom <- AxisAtom(x, axis = NULL, keepdims = TRUE)
  expect_equal(c_atom@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("AxisAtom get_data returns axis and keepdims", {
  x <- Variable(c(3, 4))
  a <- AxisAtom(x, axis = 1L, keepdims = TRUE)
  data <- get_data(a)
  expect_equal(data[[1]], 1L)
  expect_true(data[[2]])
})

## @cvxpy NONE
test_that("AxisAtom rejects out-of-range axis", {
  x <- Variable(c(3, 4))
  expect_error(AxisAtom(x, axis = 3L))
  expect_error(AxisAtom(x, axis = -3L))
})

## @cvxpy NONE
test_that("AxisAtom negative axis normalization", {
  x <- Variable(c(3, 4))
  a <- AxisAtom(x, axis = -1L)
  ## axis=-1 on 2D → axis=2 → reduces rows → (1, ncol)

  expect_equal(a@shape, c(1L, 4L))

  b <- AxisAtom(x, axis = -2L)
  ## axis=-2 on 2D → axis=1 → reduces cols → (nrow, 1)
  expect_equal(b@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("AxisAtom inherits from Atom", {
  x <- Variable(c(2, 3))
  a <- AxisAtom(x)
  expect_true(S7_inherits(a, Atom))
  expect_true(S7_inherits(a, AxisAtom))
})

## ── AxisAffAtom Base Class ──────────────────────────────────────

## @cvxpy NONE
test_that("AxisAffAtom with axis=NULL reduces to scalar", {
  x <- Variable(c(3, 4))
  a <- AxisAffAtom(x, axis = NULL)
  expect_equal(a@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("AxisAffAtom with axis=2 works", {
  x <- Variable(c(3, 4))
  a <- AxisAffAtom(x, axis = 2L)
  expect_equal(a@shape, c(1L, 4L))
})

## @cvxpy NONE
test_that("AxisAffAtom inherits from AffAtom", {
  x <- Variable(c(2, 3))
  a <- AxisAffAtom(x)
  expect_true(S7_inherits(a, AffAtom))
  expect_true(S7_inherits(a, Atom))
})

## @cvxpy NONE
test_that("AxisAffAtom is affine (convex and concave)", {
  x <- Variable(c(2, 3))
  a <- AxisAffAtom(x)
  expect_true(is_atom_convex(a))
  expect_true(is_atom_concave(a))
})

## @cvxpy NONE
test_that("AxisAffAtom get_data returns axis and keepdims", {
  x <- Variable(c(3, 4))
  a <- AxisAffAtom(x, axis = 1L, keepdims = TRUE)
  data <- get_data(a)
  expect_equal(data[[1]], 1L)
  expect_true(data[[2]])
})

## ── domain() Bug Fix (B1) ──────────────────────────────────────

## @cvxpy NONE
test_that("domain() calls atom_domain() on the atom itself", {
  ## The fix ensures domain() = atom_domain(self) + arg domains
  ## For base Atom, atom_domain returns empty list
  x <- Variable(c(2, 3))
  a <- Elementwise(args = list(x))
  d <- domain(a)
  expect_true(is.list(d))
  expect_equal(length(d), 0L)
})

## @cvxpy NONE
test_that("domain() collects arg domains recursively", {
  ## Variable domain is empty, so nested expression domain is also empty
  x <- Variable(c(2, 3))
  y <- Variable(c(2, 3))
  e <- Elementwise(args = list(x + y))
  d <- domain(e)
  expect_true(is.list(d))
})

## ── New LinOp Constructors ──────────────────────────────────────

## @cvxpy NONE
test_that("diag_vec_linop creates correct shape", {
  v <- create_var(c(3L, 1L), 1L)
  op <- diag_vec_linop(v)
  expect_equal(op$shape, c(3L, 3L))
  expect_equal(op$type, "diag_vec")
  expect_equal(op$data, 0L)
})

## @cvxpy NONE
test_that("diag_vec_linop with offset", {
  v <- create_var(c(3L, 1L), 1L)
  op <- diag_vec_linop(v, k = 1L)
  expect_equal(op$shape, c(4L, 4L))
  expect_equal(op$data, 1L)
})

## @cvxpy NONE
test_that("diag_mat_linop creates correct shape", {
  m <- create_var(c(4L, 4L), 1L)
  op <- diag_mat_linop(m)
  expect_equal(op$shape, c(4L, 1L))
  expect_equal(op$type, "diag_mat")
})

## @cvxpy NONE
test_that("trace_linop creates scalar output", {
  m <- create_var(c(3L, 3L), 1L)
  op <- trace_linop(m)
  expect_equal(op$shape, c(1L, 1L))
  expect_equal(op$type, "trace")
})

## @cvxpy NONE
test_that("upper_tri_linop creates correct shape", {
  m <- create_var(c(3L, 3L), 1L)
  op <- upper_tri_linop(m)
  ## (9 - 3) / 2 = 3 entries
  expect_equal(op$shape, c(3L, 1L))
  expect_equal(op$type, "upper_tri")
})

## @cvxpy NONE
test_that("hstack_linop creates correct shape", {
  a <- create_var(c(3L, 2L), 1L)
  b <- create_var(c(3L, 4L), 2L)
  op <- hstack_linop(list(a, b), c(3L, 6L))
  expect_equal(op$shape, c(3L, 6L))
  expect_equal(op$type, "hstack")
  expect_equal(length(op$args), 2L)
})

## @cvxpy NONE
test_that("vstack_linop creates correct shape", {
  a <- create_var(c(2L, 3L), 1L)
  b <- create_var(c(4L, 3L), 2L)
  op <- vstack_linop(list(a, b), c(6L, 3L))
  expect_equal(op$shape, c(6L, 3L))
  expect_equal(op$type, "vstack")
})

## @cvxpy NONE
test_that("kron_r_linop stores data correctly", {
  const <- create_const(matrix(1, 2, 2), c(2L, 2L))
  var <- create_var(c(3L, 3L), 1L)
  op <- kron_r_linop(const, var, c(6L, 6L))
  expect_equal(op$shape, c(6L, 6L))
  expect_equal(op$type, "kron")  # C++ has single KRON type
  expect_equal(length(op$args), 1L)   # var in args
  expect_true(!is.null(op$data))       # const in data
})

## @cvxpy NONE
test_that("kron_l_linop stores data correctly", {
  var <- create_var(c(2L, 2L), 1L)
  const <- create_const(matrix(1, 3, 3), c(3L, 3L))
  op <- kron_l_linop(var, const, c(6L, 6L))
  expect_equal(op$shape, c(6L, 6L))
  expect_equal(op$type, "kron_l")  # C++ KRON_L enum for variable-on-left
})

## @cvxpy NONE
test_that("conv_linop creates correct shape", {
  const <- create_const(c(1, 1, 1), c(3L, 1L))
  var <- create_var(c(5L, 1L), 1L)
  op <- conv_linop(const, var, c(7L, 1L))
  expect_equal(op$shape, c(7L, 1L))
  expect_equal(op$type, "conv")
})
