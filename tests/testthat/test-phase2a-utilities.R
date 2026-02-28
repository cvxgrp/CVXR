## Phase 2A: Prerequisite Utilities Tests
## Tests for sign combination, shape inference, LinOp operations, and error conditions.

# ── sign.R tests ─────────────────────────────────────────────────────

## @cvxpy NONE
test_that("sum_signs: all nonneg summands → nonneg", {
  c1 <- Constant(3)
  c2 <- Constant(5)
  result <- sum_signs(list(c1, c2))
  expect_true(result[["is_nonneg"]])
  expect_false(result[["is_nonpos"]])
})

## @cvxpy NONE
test_that("sum_signs: all nonpos summands → nonpos", {
  c1 <- Constant(-3)
  c2 <- Constant(-5)
  result <- sum_signs(list(c1, c2))
  expect_false(result[["is_nonneg"]])
  expect_true(result[["is_nonpos"]])
})

## @cvxpy NONE
test_that("sum_signs: mixed sign → unknown", {
  c1 <- Constant(3)
  c2 <- Constant(-5)
  result <- sum_signs(list(c1, c2))
  expect_false(result[["is_nonneg"]])
  expect_false(result[["is_nonpos"]])
})

## @cvxpy NONE
test_that("sum_signs: all zero → both nonneg and nonpos", {
  c1 <- Constant(0)
  c2 <- Constant(0)
  result <- sum_signs(list(c1, c2))
  expect_true(result[["is_nonneg"]])
  expect_true(result[["is_nonpos"]])
})

## @cvxpy NONE
test_that("sum_signs: single element", {
  c1 <- Constant(7)
  result <- sum_signs(list(c1))
  expect_true(result[["is_nonneg"]])
  expect_false(result[["is_nonpos"]])
})

## @cvxpy NONE
test_that("sum_signs: unknown sign variable", {
  v <- Variable(1)
  c1 <- Constant(3)
  result <- sum_signs(list(v, c1))
  ## Variable has unknown sign → sum is not nonneg

  expect_false(result[["is_nonneg"]])
  expect_false(result[["is_nonpos"]])
})

## @cvxpy NONE
test_that("sum_signs: nonneg variables", {
  v <- Variable(1, nonneg = TRUE)
  c1 <- Constant(3)
  result <- sum_signs(list(v, c1))
  expect_true(result[["is_nonneg"]])
  expect_false(result[["is_nonpos"]])
})

## @cvxpy NONE
test_that("mul_sign: nonneg * nonneg → nonneg", {
  c1 <- Constant(3)
  c2 <- Constant(5)
  result <- mul_sign(c1, c2)
  expect_true(result[["is_nonneg"]])
  expect_false(result[["is_nonpos"]])
})

## @cvxpy NONE
test_that("mul_sign: nonpos * nonpos → nonneg", {
  c1 <- Constant(-3)
  c2 <- Constant(-5)
  result <- mul_sign(c1, c2)
  expect_true(result[["is_nonneg"]])
  expect_false(result[["is_nonpos"]])
})

## @cvxpy NONE
test_that("mul_sign: nonneg * nonpos → nonpos", {
  c1 <- Constant(3)
  c2 <- Constant(-5)
  result <- mul_sign(c1, c2)
  expect_false(result[["is_nonneg"]])
  expect_true(result[["is_nonpos"]])
})

## @cvxpy NONE
test_that("mul_sign: nonpos * nonneg → nonpos", {
  c1 <- Constant(-3)
  c2 <- Constant(5)
  result <- mul_sign(c1, c2)
  expect_false(result[["is_nonneg"]])
  expect_true(result[["is_nonpos"]])
})

## @cvxpy NONE
test_that("mul_sign: zero * anything → zero (both nonneg and nonpos)", {
  c_zero <- Constant(0)
  c_pos <- Constant(5)
  result <- mul_sign(c_zero, c_pos)
  expect_true(result[["is_nonneg"]])
  expect_true(result[["is_nonpos"]])
})

## @cvxpy NONE
test_that("mul_sign: anything * zero → zero", {
  c_neg <- Constant(-5)
  c_zero <- Constant(0)
  result <- mul_sign(c_neg, c_zero)
  expect_true(result[["is_nonneg"]])
  expect_true(result[["is_nonpos"]])
})

## @cvxpy NONE
test_that("mul_sign: unknown * nonneg → unknown", {
  v <- Variable(1)
  c_pos <- Constant(5)
  result <- mul_sign(v, c_pos)
  expect_false(result[["is_nonneg"]])
  expect_false(result[["is_nonpos"]])
})

## @cvxpy NONE
test_that("mul_sign: unknown * unknown → unknown", {
  v1 <- Variable(1)
  v2 <- Variable(1)
  result <- mul_sign(v1, v2)
  expect_false(result[["is_nonneg"]])
  expect_false(result[["is_nonpos"]])
})

# ── shape.R tests ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("sum_shapes: same shapes → same result", {
  result <- sum_shapes(list(c(3L, 4L), c(3L, 4L)))
  expect_equal(result, c(3L, 4L))
})

## @cvxpy NONE
test_that("sum_shapes: broadcasting (1,n) + (m,n) → (m,n)", {
  result <- sum_shapes(list(c(1L, 3L), c(4L, 3L)))
  expect_equal(result, c(4L, 3L))
})

## @cvxpy NONE
test_that("sum_shapes: broadcasting (m,1) + (m,n) → (m,n)", {
  result <- sum_shapes(list(c(3L, 1L), c(3L, 4L)))
  expect_equal(result, c(3L, 4L))
})

## @cvxpy NONE
test_that("sum_shapes: broadcasting scalar (1,1) + (m,n) → (m,n)", {
  result <- sum_shapes(list(c(1L, 1L), c(3L, 4L)))
  expect_equal(result, c(3L, 4L))
})

## @cvxpy NONE
test_that("sum_shapes: broadcasting (1,1) + (1,1) → (1,1)", {
  result <- sum_shapes(list(c(1L, 1L), c(1L, 1L)))
  expect_equal(result, c(1L, 1L))
})

## @cvxpy NONE
test_that("sum_shapes: three shapes with broadcasting", {
  result <- sum_shapes(list(c(3L, 1L), c(1L, 4L), c(3L, 4L)))
  expect_equal(result, c(3L, 4L))
})

## @cvxpy NONE
test_that("sum_shapes: incompatible shapes → error", {
  expect_error(sum_shapes(list(c(3L, 4L), c(5L, 4L))),
               "Cannot broadcast")
})

## @cvxpy NONE
test_that("sum_shapes: single shape", {
  result <- sum_shapes(list(c(2L, 5L)))
  expect_equal(result, c(2L, 5L))
})

## @cvxpy NONE
test_that("sum_shapes: empty list → error", {
  expect_error(sum_shapes(list()), "at least one shape")
})

## @cvxpy NONE
test_that("mul_shapes: standard matmul (3,4) * (4,2) → (3,2)", {
  result <- mul_shapes(c(3L, 4L), c(4L, 2L))
  expect_equal(result, c(3L, 2L))
})

## @cvxpy NONE
test_that("mul_shapes: vector times matrix (3,1) * (1,4) → (3,4)", {
  result <- mul_shapes(c(3L, 1L), c(1L, 4L))
  expect_equal(result, c(3L, 4L))
})

## @cvxpy NONE
test_that("mul_shapes: square matmul (3,3) * (3,3) → (3,3)", {
  result <- mul_shapes(c(3L, 3L), c(3L, 3L))
  expect_equal(result, c(3L, 3L))
})

## @cvxpy NONE
test_that("mul_shapes: (1,1) * (1,n) → (1,n)", {
  result <- mul_shapes(c(1L, 1L), c(1L, 5L))
  expect_equal(result, c(1L, 5L))
})

## @cvxpy NONE
test_that("mul_shapes: incompatible inner dimensions → error", {
  expect_error(mul_shapes(c(3L, 4L), c(5L, 2L)),
               "Incompatible dimensions")
})

## @cvxpy NONE
test_that("mul_shapes_promote: returns list with correct components", {
  result <- mul_shapes_promote(c(3L, 4L), c(4L, 2L))
  expect_type(result, "list")
  expect_equal(result$lh_shape, c(3L, 4L))
  expect_equal(result$rh_shape, c(4L, 2L))
  expect_equal(result$shape, c(3L, 2L))
})

# ── error.R tests ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("DCPError creates condition of correct class", {
  err <- DCPError("test message")
  expect_s3_class(err, "DCPError")
  expect_s3_class(err, "error")
  expect_s3_class(err, "condition")
  expect_equal(err$message, "test message")
})

## @cvxpy NONE
test_that("DCPError can be caught with tryCatch", {
  caught <- FALSE
  tryCatch(
    stop(DCPError("dcp violation")),
    DCPError = function(e) {
      caught <<- TRUE
      expect_equal(e$message, "dcp violation")
    }
  )
  expect_true(caught)
})

## @cvxpy NONE
test_that("SolverError creates condition of correct class", {
  err <- SolverError("solver failed")
  expect_s3_class(err, "SolverError")
  expect_s3_class(err, "error")
  expect_equal(err$message, "solver failed")
})

## @cvxpy NONE
test_that("SolverError can be caught with tryCatch", {
  caught <- FALSE
  tryCatch(
    stop(SolverError("solver failed")),
    SolverError = function(e) {
      caught <<- TRUE
    }
  )
  expect_true(caught)
})

# ── LinOp operation function tests ───────────────────────────────────

## @cvxpy NONE
test_that("is_scalar_linop: scalar constant is scalar", {
  op <- create_const(5, c(1L, 1L))
  expect_true(CVXR:::is_scalar_linop(op))
})

## @cvxpy NONE
test_that("is_scalar_linop: matrix is not scalar", {
  op <- create_const(matrix(1:6, 2, 3), c(2L, 3L))
  expect_false(CVXR:::is_scalar_linop(op))
})

## @cvxpy NONE
test_that("is_const_linop: scalar const is constant", {
  op <- create_const(5, c(1L, 1L))
  expect_true(CVXR:::is_const_linop(op))
})

## @cvxpy NONE
test_that("is_const_linop: variable is not constant", {
  op <- create_var(c(3L, 1L), 1L)
  expect_false(CVXR:::is_const_linop(op))
})

## @cvxpy NONE
test_that("is_const_linop: param is constant", {
  op <- create_param(c(3L, 1L), 1L)
  expect_true(CVXR:::is_const_linop(op))
})

## @cvxpy NONE
test_that("sum_expr_linop: creates SUM LinOp with correct type and shape", {
  op1 <- create_var(c(3L, 1L), 1L)
  op2 <- create_var(c(3L, 1L), 2L)
  result <- CVXR:::sum_expr_linop(list(op1, op2))
  expect_equal(result$type, "sum")
  expect_equal(result$shape, c(3L, 1L))
  expect_length(result$args, 2L)
  expect_null(result$data)
})

## @cvxpy NONE
test_that("neg_expr_linop: creates NEG LinOp", {
  op <- create_var(c(3L, 1L), 1L)
  result <- CVXR:::neg_expr_linop(op)
  expect_equal(result$type, "neg")
  expect_equal(result$shape, c(3L, 1L))
  expect_length(result$args, 1L)
  expect_identical(result$args[[1L]], op)
  expect_null(result$data)
})

## @cvxpy NONE
test_that("mul_expr_linop: creates MUL LinOp with const as data", {
  const_op <- create_const(matrix(1:6, 2, 3), c(2L, 3L))
  var_op <- create_var(c(3L, 1L), 1L)
  result <- CVXR:::mul_expr_linop(const_op, var_op, c(2L, 1L))
  expect_equal(result$type, "mul_expr")
  expect_equal(result$shape, c(2L, 1L))
  expect_length(result$args, 1L)
  expect_identical(result$args[[1L]], var_op)
  expect_identical(result$data, const_op)
})

## @cvxpy NONE
test_that("rmul_expr_linop: creates RMUL LinOp with const as data", {
  var_op <- create_var(c(2L, 3L), 1L)
  const_op <- create_const(matrix(1:12, 3, 4), c(3L, 4L))
  result <- CVXR:::rmul_expr_linop(var_op, const_op, c(2L, 4L))
  expect_equal(result$type, "rmul_expr")
  expect_equal(result$shape, c(2L, 4L))
  expect_length(result$args, 1L)
  expect_identical(result$args[[1L]], var_op)
  expect_identical(result$data, const_op)
})

## @cvxpy NONE
test_that("multiply_linop: creates MUL_ELEM LinOp with broadcast shape", {
  const_op <- create_const(matrix(1:3, 3, 1), c(3L, 1L))
  var_op <- create_var(c(3L, 1L), 1L)
  result <- CVXR:::multiply_linop(const_op, var_op)
  expect_equal(result$type, "mul_elem")
  expect_equal(result$shape, c(3L, 1L))
  expect_length(result$args, 1L)
  expect_identical(result$args[[1L]], var_op)
  expect_identical(result$data, const_op)
})

## @cvxpy NONE
test_that("div_expr_linop: creates DIV LinOp", {
  var_op <- create_var(c(3L, 1L), 1L)
  const_op <- create_const(2, c(1L, 1L))
  result <- CVXR:::div_expr_linop(var_op, const_op)
  expect_equal(result$type, "div")
  expect_equal(result$shape, c(3L, 1L))
  expect_length(result$args, 1L)
  expect_identical(result$args[[1L]], var_op)
  expect_identical(result$data, const_op)
})

## @cvxpy NONE
test_that("promote_linop: creates PROMOTE LinOp", {
  scalar_op <- create_const(5, c(1L, 1L))
  result <- CVXR:::promote_linop(scalar_op, c(3L, 4L))
  expect_equal(result$type, "promote")
  expect_equal(result$shape, c(3L, 4L))
  expect_length(result$args, 1L)
  expect_identical(result$args[[1L]], scalar_op)
  expect_null(result$data)
})

## @cvxpy NONE
test_that("transpose_linop: creates TRANSPOSE LinOp with reversed shape", {
  op <- create_var(c(3L, 4L), 1L)
  result <- CVXR:::transpose_linop(op)
  expect_equal(result$type, "transpose")
  expect_equal(result$shape, c(4L, 3L))
  expect_length(result$args, 1L)
  expect_identical(result$args[[1L]], op)
})

## @cvxpy NONE
test_that("index_linop: creates INDEX LinOp with keys as data", {
  op <- create_var(c(3L, 4L), 1L)
  keys <- list(1:2, 2:3)
  result <- CVXR:::index_linop(op, c(2L, 2L), keys)
  expect_equal(result$type, "index")
  expect_equal(result$shape, c(2L, 2L))
  expect_length(result$args, 1L)
  ## index_linop appends "key" marker for C++ bridge (set_slice_data)
  expect_identical(result$data, c(keys, list("key")))
})

## @cvxpy NONE
test_that("sum_entries_linop: creates SUM_ENTRIES LinOp", {
  op <- create_var(c(3L, 4L), 1L)
  result <- CVXR:::sum_entries_linop(op, c(1L, 1L))
  expect_equal(result$type, "sum_entries")
  expect_equal(result$shape, c(1L, 1L))
  expect_length(result$args, 1L)
})

## @cvxpy NONE
test_that("reshape_linop: creates RESHAPE LinOp", {
  op <- create_var(c(3L, 4L), 1L)
  result <- CVXR:::reshape_linop(op, c(12L, 1L))
  expect_equal(result$type, "reshape_expr")
  expect_equal(result$shape, c(12L, 1L))
  expect_length(result$args, 1L)
})

# ── LinOp type constants existence check ─────────────────────────────

## @cvxpy NONE
test_that("All LinOp type constants are defined", {
  expect_equal(CVXR:::LINOP_PROMOTE, "promote")
  expect_equal(CVXR:::LINOP_MUL, "mul_expr")
  expect_equal(CVXR:::LINOP_RMUL, "rmul_expr")
  expect_equal(CVXR:::LINOP_MUL_ELEM, "mul_elem")
  expect_equal(CVXR:::LINOP_DIV, "div")
  expect_equal(CVXR:::LINOP_SUM, "sum")
  expect_equal(CVXR:::LINOP_NEG, "neg")
  expect_equal(CVXR:::LINOP_INDEX, "index")
  expect_equal(CVXR:::LINOP_TRANSPOSE, "transpose")
  expect_equal(CVXR:::LINOP_SUM_ENTRIES, "sum_entries")
  expect_equal(CVXR:::LINOP_RESHAPE, "reshape_expr")
})

## @cvxpy NONE
test_that("LinOp type strings match CVXPY lin_op.py exactly", {
  ## CRITICAL: These strings must match CVXPY exactly because
  ## the C++ canonicalization backend uses them.
  expect_equal(CVXR:::LINOP_VARIABLE, "variable")
  expect_equal(CVXR:::LINOP_PARAM, "param")
  expect_equal(CVXR:::LINOP_SCALAR_CONST, "scalar_const")
  expect_equal(CVXR:::LINOP_DENSE_CONST, "dense_const")
  expect_equal(CVXR:::LINOP_SPARSE_CONST, "sparse_const")
  expect_equal(CVXR:::LINOP_NO_OP, "no_op")
})
