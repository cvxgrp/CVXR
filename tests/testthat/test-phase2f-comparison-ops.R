## Tests for Phase 2F: Comparison Operator Dispatch (==, <=, >=)

# ── == creates Equality ───────────────────────────────────────────────

## @cvxpy NONE
test_that("Variable == Variable creates Equality", {
  x <- Variable(3, name = "x")
  y <- Variable(3, name = "y")
  constr <- x == y
  expect_s3_class(constr, "CVXR::Equality")
  expect_equal(constr@shape, c(3L, 1L))
  expect_length(constr@args, 2L)
})

## @cvxpy NONE
test_that("Variable == numeric creates Equality", {
  x <- Variable(3)
  constr <- x == 0
  expect_s3_class(constr, "CVXR::Equality")
  expect_equal(constr@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("numeric == Variable creates Equality (reverse dispatch)", {
  x <- Variable(3)
  constr <- 0 == x
  expect_s3_class(constr, "CVXR::Equality")
  expect_equal(constr@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Variable == matrix creates Equality", {
  x <- Variable(c(2, 3))
  m <- matrix(1:6, 2, 3)
  constr <- x == m
  expect_s3_class(constr, "CVXR::Equality")
  expect_equal(constr@shape, c(2L, 3L))
})

## @cvxpy NONE
test_that("== Equality is_dcp when affine", {
  x <- Variable(3)
  constr <- x == 0
  expect_true(is_dcp(constr))
})

# ── <= creates Inequality ─────────────────────────────────────────────

## @cvxpy NONE
test_that("Variable <= Variable creates Inequality", {
  x <- Variable(3, name = "x")
  y <- Variable(3, name = "y")
  constr <- x <= y
  expect_s3_class(constr, "CVXR::Inequality")
  expect_equal(constr@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Variable <= numeric creates Inequality", {
  x <- Variable(3)
  constr <- x <= 10
  expect_s3_class(constr, "CVXR::Inequality")
  expect_equal(constr@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("numeric <= Variable creates Inequality (reverse)", {
  x <- Variable(3)
  constr <- 0 <= x
  expect_s3_class(constr, "CVXR::Inequality")
  ## 0 <= x  ⟹  Inequality(x_rhs=Constant(0), lhs=x)
  ## Actually: 0 <= x  ⟹  >= handler ⟹  Inequality(e2=Constant(0), e1=x)
  ## Wait: .Generic when 0 <= Variable is "<=", with e1=0, e2=Variable
  ## But chooseOpsMethod makes us win, so e1=0, e2=Variable, .Generic="<="
  ## Actually R will try e2's method when e1 doesn't match, with .Generic="<="
  ## and reverse=TRUE? No — R's Ops dispatch: if right operand wins,

  ## .Generic stays the same but e1/e2 are NOT swapped.
  ## So: 0 <= x → .Generic="<=", e1=Constant(0), e2=x → Inequality(Constant(0), x)
  ## That means: Constant(0) - x <= 0, i.e. x >= 0. Correct!
  expect_equal(constr@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("<= is_dcp when lhs - rhs is convex", {
  x <- Variable(3)
  constr <- x <= 10
  expect_true(is_dcp(constr))
})

# ── >= creates Inequality (reversed) ─────────────────────────────────

## @cvxpy NONE
test_that("Variable >= Variable creates Inequality", {
  x <- Variable(3, name = "x")
  y <- Variable(3, name = "y")
  constr <- x >= y
  expect_s3_class(constr, "CVXR::Inequality")
  expect_equal(constr@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Variable >= numeric creates Inequality", {
  x <- Variable(3)
  constr <- x >= 0
  expect_s3_class(constr, "CVXR::Inequality")
  ## x >= 0  ⟹  .Generic=">=", Inequality(e2=x, e1=Constant(0))
  ## = Inequality(Constant(0), x) = 0 - x <= 0 = x >= 0. Correct!
  expect_equal(constr@shape, c(3L, 1L))
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("numeric >= Variable creates Inequality (reverse)", {
  x <- Variable(3)
  constr <- 10 >= x
  expect_s3_class(constr, "CVXR::Inequality")
  expect_equal(constr@shape, c(3L, 1L))
})

# ── Value semantics through constraints ───────────────────────────────

## @cvxpy NONE
test_that("== constraint value works with reference semantics", {
  x <- Variable(c(2, 1))
  z <- Variable(c(2, 1))
  constr <- x == z
  val <- matrix(c(2, 2), 2, 1)
  value(x) <- val
  value(z) <- val
  expect_true(value(constr))
})

## @cvxpy NONE
test_that("<= constraint value works with reference semantics", {
  x <- Variable(c(2, 1))
  constr <- x <= 10
  value(x) <- matrix(c(3, 5), 2, 1)
  expect_true(value(constr))
})

## @cvxpy NONE
test_that(">= constraint value works with reference semantics", {
  x <- Variable(c(2, 1))
  constr <- x >= 0
  value(x) <- matrix(c(3, 5), 2, 1)
  expect_true(value(constr))
})

# ── DCP properties ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("<= with affine expression is DCP", {
  x <- Variable(3)
  y <- Variable(3)
  constr <- (x + y) <= 10
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that(">= with affine expression is DCP", {
  x <- Variable(3)
  y <- Variable(3)
  constr <- (x + y) >= 0
  expect_true(is_dcp(constr))
})

## @cvxpy NONE
test_that("== with affine expression is DCP", {
  x <- Variable(3)
  y <- Variable(3)
  constr <- (x - y) == 0
  expect_true(is_dcp(constr))
})

# ── Constraint list from multiple comparisons ─────────────────────────

## @cvxpy NONE
test_that("can build constraint list with operators", {
  x <- Variable(3)
  constraints <- list(x >= 0, x <= 10)
  expect_length(constraints, 2L)
  expect_s3_class(constraints[[1]], "CVXR::Inequality")
  expect_s3_class(constraints[[2]], "CVXR::Inequality")
  expect_true(all(vapply(constraints, is_dcp, logical(1))))
})

# ── Matrix operands ───────────────────────────────────────────────────

## @cvxpy NONE
test_that("Variable == matrix works", {
  x <- Variable(c(2, 3))
  m <- matrix(1:6, 2, 3)
  constr <- x == m
  expect_s3_class(constr, "CVXR::Equality")
  expect_equal(constr@shape, c(2L, 3L))
})

## @cvxpy NONE
test_that("Variable <= matrix works", {
  x <- Variable(c(2, 3))
  m <- matrix(1:6, 2, 3)
  constr <- x <= m
  expect_s3_class(constr, "CVXR::Inequality")
  expect_equal(constr@shape, c(2L, 3L))
})

## @cvxpy NONE
test_that("Variable >= matrix works", {
  x <- Variable(c(2, 3))
  m <- matrix(0, 2, 3)
  constr <- x >= m
  expect_s3_class(constr, "CVXR::Inequality")
  expect_equal(constr@shape, c(2L, 3L))
})

# ── Scalar comparisons ───────────────────────────────────────────────

## @cvxpy NONE
test_that("scalar Variable == scalar numeric", {
  x <- Variable(1)
  constr <- x == 5
  expect_s3_class(constr, "CVXR::Equality")
  expect_equal(constr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("scalar Variable <= scalar numeric", {
  x <- Variable(1)
  constr <- x <= 5
  expect_s3_class(constr, "CVXR::Inequality")
  expect_equal(constr@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("scalar Variable >= scalar numeric", {
  x <- Variable(1)
  constr <- x >= 0
  expect_s3_class(constr, "CVXR::Inequality")
  expect_equal(constr@shape, c(1L, 1L))
})

# ── Expression operands (not just leaf) ──────────────────────────────

## @cvxpy NONE
test_that("expression <= expression works", {
  x <- Variable(3)
  y <- Variable(3)
  constr <- (x + y) <= (2 * x)
  expect_s3_class(constr, "CVXR::Inequality")
  expect_equal(constr@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("expression == expression works", {
  x <- Variable(3)
  y <- Variable(3)
  constr <- (x + y) == (2 * x)
  expect_s3_class(constr, "CVXR::Equality")
  expect_equal(constr@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("expression >= expression works", {
  x <- Variable(3)
  y <- Variable(3)
  constr <- (2 * x) >= (x + y)
  expect_s3_class(constr, "CVXR::Inequality")
  expect_equal(constr@shape, c(3L, 1L))
})
