## Tests for Phase 2G: Objective Classes (Minimize, Maximize)
## Translated from CVXPY: tests/test_objectives.py

# ── Minimize construction ─────────────────────────────────────────────

## @cvxpy NONE
test_that("Minimize accepts scalar variable", {
  x <- Variable(1, name = "x")
  obj <- Minimize(x)
  expect_s3_class(obj, "CVXR::Minimize")
  expect_length(obj@args, 1L)
})

## @cvxpy NONE
test_that("Minimize accepts scalar expression", {
  x <- Variable(1, name = "x")
  z <- Variable(1, name = "z")
  obj <- Minimize(x + z)
  expect_s3_class(obj, "CVXR::Minimize")
})

## @cvxpy test_objectives.py::TestObjectives::test_minimize
test_that("Minimize rejects non-scalar", {
  y <- Variable(3, name = "y")
  expect_error(Minimize(y), "must resolve to a scalar")
})

## @cvxpy NONE
test_that("Minimize accepts scalar constant", {
  obj <- Minimize(Constant(5))
  expect_s3_class(obj, "CVXR::Minimize")
})

## @cvxpy NONE
test_that("Minimize accepts numeric (auto-promoted)", {
  obj <- Minimize(5)
  expect_s3_class(obj, "CVXR::Minimize")
})

# ── Maximize construction ─────────────────────────────────────────────

## @cvxpy NONE
test_that("Maximize accepts scalar variable", {
  x <- Variable(1, name = "x")
  obj <- Maximize(x)
  expect_s3_class(obj, "CVXR::Maximize")
  expect_length(obj@args, 1L)
})

## @cvxpy NONE
test_that("Maximize accepts scalar expression", {
  x <- Variable(1, name = "x")
  z <- Variable(1, name = "z")
  obj <- Maximize(x + z)
  expect_s3_class(obj, "CVXR::Maximize")
})

## @cvxpy test_objectives.py::TestObjectives::test_maximize
test_that("Maximize rejects non-scalar", {
  y <- Variable(3, name = "y")
  expect_error(Maximize(y), "must resolve to a scalar")
})

# ── DCP checks ────────────────────────────────────────────────────────
## CVXPY: test_objectives.py::test_is_dcp

## @cvxpy test_objectives.py::TestObjectives::test_is_dcp
test_that("Minimize of affine is DCP", {
  x <- Variable(1, name = "x")
  obj <- Minimize(x)
  expect_true(is_dcp(obj))
})

## @cvxpy test_objectives.py::TestObjectives::test_is_dcp
test_that("Minimize of constant is DCP", {
  obj <- Minimize(Constant(5))
  expect_true(is_dcp(obj))
})

## @cvxpy test_objectives.py::TestObjectives::test_is_dcp
test_that("Minimize of convex expression is DCP", {
  ## x + z is affine (hence convex)
  x <- Variable(1, name = "x")
  z <- Variable(1, name = "z")
  obj <- Minimize(x + z)
  expect_true(is_dcp(obj))
})

## @cvxpy test_objectives.py::TestObjectives::test_is_dcp
test_that("Maximize of affine is DCP", {
  x <- Variable(1, name = "x")
  obj <- Maximize(x)
  expect_true(is_dcp(obj))
})

## @cvxpy test_objectives.py::TestObjectives::test_is_dcp
test_that("Maximize of constant is DCP", {
  obj <- Maximize(Constant(5))
  expect_true(is_dcp(obj))
})

## @cvxpy test_objectives.py::TestObjectives::test_is_dcp
test_that("Maximize of concave expression is DCP (affine is concave)", {
  x <- Variable(1, name = "x")
  z <- Variable(1, name = "z")
  obj <- Maximize(x + z)
  expect_true(is_dcp(obj))
})

# ── Value ─────────────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Minimize value extracts scalar", {
  obj <- Minimize(Constant(5))
  expect_equal(value(obj), 5)
})

## @cvxpy NONE
test_that("Maximize value extracts scalar", {
  obj <- Maximize(Constant(3))
  expect_equal(value(obj), 3)
})

## @cvxpy NONE
test_that("Minimize value is NULL when variable has no value", {
  x <- Variable(1)
  obj <- Minimize(x)
  expect_null(value(obj))
})

## @cvxpy NONE
test_that("Maximize value is NULL when variable has no value", {
  x <- Variable(1)
  obj <- Maximize(x)
  expect_null(value(obj))
})

## @cvxpy NONE
test_that("Minimize value with assigned variable", {
  x <- Variable(1)
  obj <- Minimize(x)
  value(x) <- 7
  expect_equal(value(obj), 7)
})

## @cvxpy NONE
test_that("Maximize value with assigned variable", {
  x <- Variable(1)
  obj <- Maximize(x)
  value(x) <- 3
  expect_equal(value(obj), 3)
})

# ── Print / name ──────────────────────────────────────────────────────
## CVXPY: test_objectives.py::test_str

## @cvxpy test_objectives.py::TestObjectives::test_str
test_that("Minimize print includes 'minimize'", {
  x <- Variable(1, name = "x")
  obj <- Minimize(x)
  expect_output(print(obj), "minimize")
})

## @cvxpy test_objectives.py::TestObjectives::test_str
test_that("Maximize print includes 'maximize'", {
  x <- Variable(1, name = "x")
  obj <- Maximize(x)
  expect_output(print(obj), "maximize")
})

# ── Canonicalize ──────────────────────────────────────────────────────
## CVXPY: test_objectives.py::test_minimize (canonical_form)

## @cvxpy test_objectives.py::TestObjectives::test_minimize
test_that("Minimize canonicalize for affine has no constraints", {
  x <- Variable(1, name = "x")
  z <- Variable(1, name = "z")
  exp <- x + z
  obj <- Minimize(exp)
  result <- canonical_form(obj)
  expect_length(result, 2L)
  ## For affine objectives, there should be no constraints
  expect_length(result[[2L]], 0L)
})

## @cvxpy test_objectives.py::TestObjectives::test_maximize
test_that("Maximize canonicalize for affine has no constraints", {
  x <- Variable(1, name = "x")
  z <- Variable(1, name = "z")
  exp <- x + z
  obj <- Maximize(exp)
  result <- canonical_form(obj)
  expect_length(result, 2L)
  expect_length(result[[2L]], 0L)
})

# ── Copy ──────────────────────────────────────────────────────────────
## CVXPY: test_objectives.py::test_minimize (copy)

## @cvxpy test_objectives.py::TestObjectives::test_minimize
test_that("Minimize copy creates new object", {
  x <- Variable(1, name = "x")
  z <- Variable(1, name = "z")
  obj <- Minimize(x + z)
  copy_obj <- expr_copy(obj)
  expect_s3_class(copy_obj, "CVXR::Minimize")
  ## New object with same class
  expect_false(identical(copy_obj@id, obj@id))
})

## @cvxpy test_objectives.py::TestObjectives::test_maximize
test_that("Maximize copy creates new object", {
  x <- Variable(1, name = "x")
  z <- Variable(1, name = "z")
  obj <- Maximize(x + z)
  copy_obj <- expr_copy(obj)
  expect_s3_class(copy_obj, "CVXR::Maximize")
  expect_false(identical(copy_obj@id, obj@id))
})

# ── Variables / Parameters / Constants ────────────────────────────────

## @cvxpy NONE
test_that("Minimize variables returns contained variables", {
  x <- Variable(1, name = "x")
  z <- Variable(1, name = "z")
  obj <- Minimize(x + z)
  vars <- variables(obj)
  expect_length(vars, 2L)
})

## @cvxpy NONE
test_that("Maximize parameters returns contained parameters", {
  p <- Parameter(1, name = "p")
  x <- Variable(1, name = "x")
  obj <- Maximize(x + p)
  params <- parameters(obj)
  expect_length(params, 1L)
})

## @cvxpy NONE
test_that("Minimize constants returns contained constants", {
  x <- Variable(1, name = "x")
  obj <- Minimize(x + 5)
  consts <- constants(obj)
  expect_length(consts, 1L)
})

# ── scalar_value utility ──────────────────────────────────────────────

## @cvxpy NONE
test_that("scalar_value extracts from 1x1 matrix", {
  expect_equal(scalar_value(matrix(5, 1, 1)), 5)
})

## @cvxpy NONE
test_that("scalar_value handles plain numeric", {
  expect_equal(scalar_value(42), 42)
})

## @cvxpy NONE
test_that("scalar_value handles NULL", {
  expect_null(scalar_value(NULL))
})

## @cvxpy NONE
test_that("scalar_value extracts from 1x1 sparse matrix", {
  m <- Matrix::sparseMatrix(i = 1, j = 1, x = 3, dims = c(1, 1))
  expect_equal(scalar_value(m), 3)
})

# ── Objective arithmetic ────────────────────────────────────────────

## @cvxpy test_objectives.py::TestObjectives::test_add_problems
test_that("Objective arithmetic: add, subtract, multiply (CVXPY parity)", {
  x <- Variable(1, name = "x")
  expr1 <- power(x, 2)
  expr2 <- power(x, -1)
  alpha <- 2

  ## Minimize + Minimize is DCP
  sum_min <- CVXR:::.add_objectives(Minimize(expr1), Minimize(expr2))
  expect_true(is_dcp(sum_min))
  expect_s3_class(sum_min, "CVXR::Minimize")

  ## Maximize + Maximize is DCP (concave + concave = concave)
  sum_max <- CVXR:::.add_objectives(Maximize(-expr1), Maximize(-expr2))
  expect_true(is_dcp(sum_max))
  expect_s3_class(sum_max, "CVXR::Maximize")

  ## Minimize + Maximize raises DCP error
  expect_error(
    CVXR:::.add_objectives(Minimize(expr1), Maximize(-expr2)),
    "DCP"
  )

  ## Minimize - Maximize is DCP (becomes Minimize + Minimize)
  diff_obj <- CVXR:::.sub_objectives(Minimize(expr1), Maximize(-expr2))
  expect_true(is_dcp(diff_obj))
  expect_s3_class(diff_obj, "CVXR::Minimize")

  ## alpha * Minimize is DCP (positive scalar preserves direction)
  scaled_min <- CVXR:::.mul_objective(Minimize(expr1), alpha)
  expect_true(is_dcp(scaled_min))

  ## alpha * Maximize is DCP (positive scalar preserves direction)
  scaled_max <- CVXR:::.mul_objective(Maximize(-expr1), alpha)
  expect_true(is_dcp(scaled_max))

  ## -alpha * Maximize is DCP (negative scalar flips to Minimize)
  neg_scaled <- CVXR:::.mul_objective(Maximize(-expr1), -alpha)
  expect_true(is_dcp(neg_scaled))
  expect_s3_class(neg_scaled, "CVXR::Minimize")
})
