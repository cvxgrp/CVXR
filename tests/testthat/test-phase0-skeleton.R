## Phase 0: Skeleton tests
## Verify package loads and core infrastructure works

## @cvxpy NONE
test_that("package loads without error", {
  # S7 methods registered — .onLoad is in the package namespace
  expect_true(is.function(CVXR:::.onLoad))
})

## @cvxpy NONE
test_that("expression ID generator works", {
  CVXR:::reset_expr_id()
  id1 <- CVXR:::next_expr_id()
  id2 <- CVXR:::next_expr_id()
  expect_equal(id1, 1L)
  expect_equal(id2, 2L)
  expect_true(id2 > id1)
  CVXR:::reset_expr_id()
})

## @cvxpy NONE
test_that("S7 generics are defined", {
  is_s7_generic <- function(x) inherits(x, "S7_generic")
  # Expression interface
  expect_true(is_s7_generic(value))
  expect_true(is_s7_generic(is_convex))
  expect_true(is_s7_generic(is_concave))
  expect_true(is_s7_generic(is_dcp))
  expect_true(is_s7_generic(canonicalize))
  expect_true(is_s7_generic(expr_sign))
  # Sign/attribute queries
  expect_true(is_s7_generic(is_nonneg))
  expect_true(is_s7_generic(is_nonpos))
  expect_true(is_s7_generic(is_zero))
  expect_true(is_s7_generic(is_psd))
  expect_true(is_s7_generic(is_nsd))
  expect_true(is_s7_generic(is_symmetric))
  # Atom hooks
  expect_true(is_s7_generic(shape_from_args))
  expect_true(is_s7_generic(sign_from_args))
  expect_true(is_s7_generic(is_atom_convex))
  expect_true(is_s7_generic(is_atom_concave))
  expect_true(is_s7_generic(is_incr))
  expect_true(is_s7_generic(is_decr))
  expect_true(is_s7_generic(graph_implementation))
  expect_true(is_s7_generic(numeric_value))
})

## @cvxpy NONE
test_that("generics with extra formals have correct signatures", {
  # B1/H7/H8 fix: verify multi-param generics
  expect_true("idx" %in% names(formals(is_incr)))
  expect_true("idx" %in% names(formals(is_decr)))
  expect_true("arg_objs" %in% names(formals(graph_implementation)))
  expect_true("shape" %in% names(formals(graph_implementation)))
  expect_true("data" %in% names(formals(graph_implementation)))
  expect_true("values" %in% names(formals(numeric_value)))
})

## @cvxpy NONE
test_that("cache helpers work with actual S7 object", {
  # M6 fix: test cache_get/set/has/clear on an actual S7 object
  CacheTest <- S7::new_class("CacheTest", package = "CVXR",
    properties = list(.cache = S7::new_property(class = S7::class_any,
      default = quote(new.env(parent = emptyenv())))))
  obj <- CacheTest()

  expect_false(CVXR:::cache_has(obj, "foo"))
  CVXR:::cache_set(obj, "foo", 42)
  expect_true(CVXR:::cache_has(obj, "foo"))
  expect_equal(CVXR:::cache_get(obj, "foo"), 42)
  expect_true(CVXR:::cache_miss(CVXR:::cache_get(obj, "bar")))
  CVXR:::cache_clear(obj)
  expect_false(CVXR:::cache_has(obj, "foo"))
})

## @cvxpy NONE
test_that("cache helpers work", {
  # Create a mock object with a .cache environment
  cache <- new.env(parent = emptyenv())

  expect_false(exists("foo", envir = cache, inherits = FALSE))
  assign("foo", 42, envir = cache)
  expect_true(exists("foo", envir = cache, inherits = FALSE))
  expect_equal(get("foo", envir = cache, inherits = FALSE), 42)
  rm("foo", envir = cache)
  expect_false(exists("foo", envir = cache, inherits = FALSE))
})

## @cvxpy NONE
test_that("settings constants are defined", {
  expect_equal(SCS_SOLVER, "SCS")
  expect_equal(OSQP_SOLVER, "OSQP")
  expect_equal(CLARABEL_SOLVER, "CLARABEL")
  expect_equal(AFFINE, "AFFINE")
  expect_equal(CONVEX, "CONVEX")
  expect_equal(CONCAVE, "CONCAVE")
  expect_equal(ZERO_SIGN, "ZERO")
  expect_equal(NONNEG_SIGN, "NONNEGATIVE")
  expect_equal(NONPOS_SIGN, "NONPOSITIVE")
  expect_equal(UNKNOWN_SIGN, "UNKNOWN")
})

## @cvxpy NONE
test_that("validate_shape works correctly", {
  expect_equal(CVXR:::validate_shape(NULL), c(1L, 1L))
  expect_equal(CVXR:::validate_shape(5), c(5L, 1L))
  expect_equal(CVXR:::validate_shape(c(3, 4)), c(3L, 4L))
  expect_error(CVXR:::validate_shape(c(1, 2, 3)), "length 1 or 2")
  expect_error(CVXR:::validate_shape(c(-1, 2)), "positive")
  expect_error(CVXR:::validate_shape(c(0, 2)), "positive")
})

## @cvxpy NONE
test_that("is_scalar_shape works correctly", {
  expect_true(CVXR:::is_scalar_shape(c(1L, 1L)))
  expect_false(CVXR:::is_scalar_shape(c(3L, 1L)))
  expect_false(CVXR:::is_scalar_shape(c(1L, 3L)))
})
