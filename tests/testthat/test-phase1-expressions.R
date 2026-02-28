## Phase 1: Core Expression System tests
## Tests for Variable, Constant, Parameter, Expression, Leaf, Canonical

# ── Variable tests ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Variable construction with scalar shape", {
  CVXR:::reset_expr_id()
  x <- Variable(3)
  expect_equal(x@shape, c(3L, 1L))
  expect_equal(expr_name(x), "var1")
})

## @cvxpy NONE
test_that("Variable construction with 2D shape", {
  x <- Variable(c(2, 3))
  expect_equal(x@shape, c(2L, 3L))
})

## @cvxpy NONE
test_that("Variable construction with custom name", {
  x <- Variable(3, name = "myvar")
  expect_equal(expr_name(x), "myvar")
})

## @cvxpy NONE
test_that("Variable rejects non-string name", {
  expect_error(Variable(3, name = 42), "must be a string")
})

## @cvxpy NONE
test_that("Variable unique IDs", {
  CVXR:::reset_expr_id()
  x1 <- Variable(2)
  x2 <- Variable(3)
  expect_true(x1@id != x2@id)
})

## @cvxpy NONE
test_that("Variable custom ID", {
  x <- Variable(3, var_id = 999L)
  expect_equal(x@id, 999L)
  expect_equal(expr_name(x), "var999")
})

## @cvxpy NONE
test_that("Variable is not constant", {
  x <- Variable(3)
  expect_false(is_constant(x))
})

## @cvxpy NONE
test_that("Variable is affine (convex and concave)", {
  x <- Variable(3)
  expect_true(is_convex(x))
  expect_true(is_concave(x))
  expect_true(is_affine(x))
  expect_true(is_dcp(x))
})

## @cvxpy NONE
test_that("Variable sign: default is unknown", {
  x <- Variable(3)
  expect_false(is_nonneg(x))
  expect_false(is_nonpos(x))
  expect_false(is_zero(x))
  expect_equal(expr_sign(x), UNKNOWN_SIGN)
})

## @cvxpy NONE
test_that("Variable nonneg attribute", {
  x <- Variable(3, nonneg = TRUE)
  expect_true(is_nonneg(x))
  expect_false(is_nonpos(x))
})

## @cvxpy NONE
test_that("Variable nonpos attribute", {
  x <- Variable(3, nonpos = TRUE)
  expect_false(is_nonneg(x))
  expect_true(is_nonpos(x))
})

## @cvxpy NONE
test_that("Variable PSD attribute requires square", {
  expect_error(Variable(c(2, 3), PSD = TRUE), "square matrix")
  x <- Variable(c(3, 3), PSD = TRUE)
  expect_true(is_psd(x))
  expect_true(is_symmetric(x))
})

## @cvxpy NONE
test_that("Variable boolean attribute", {
  x <- Variable(3, boolean = TRUE)
  expect_true(is_nonneg(x))
})

## @cvxpy NONE
test_that("Variable symmetric attribute requires square", {
  expect_error(Variable(c(2, 3), symmetric = TRUE), "square matrix")
  x <- Variable(c(3, 3), symmetric = TRUE)
  expect_true(is_symmetric(x))
})

## @cvxpy NONE
test_that("Variable variables() returns self", {
  x <- Variable(3)
  vars <- variables(x)
  expect_length(vars, 1)
  expect_identical(vars[[1]]@id, x@id)
})

## @cvxpy NONE
test_that("Variable parameters() returns empty", {
  x <- Variable(3)
  expect_length(parameters(x), 0)
})

## @cvxpy NONE
test_that("Variable constants() returns empty", {
  x <- Variable(3)
  expect_length(constants(x), 0)
})

## @cvxpy NONE
test_that("Variable atoms() returns empty", {
  x <- Variable(3)
  expect_length(atoms(x), 0)
})

## @cvxpy NONE
test_that("Variable grad returns identity", {
  x <- Variable(3)
  g <- grad(x)
  expect_length(g, 1)
  key <- as.character(x@id)
  expect_true(key %in% names(g))
  mat <- g[[key]]
  expect_equal(dim(mat), c(3L, 3L))
  expect_equal(as.matrix(mat), diag(3))
})

## @cvxpy NONE
test_that("Variable canonicalize", {
  x <- Variable(c(2, 3))
  canon <- canonicalize(x)
  expect_length(canon, 2)
  linop <- canon[[1]]
  expect_equal(linop$type, "variable")
  expect_equal(linop$shape, c(2L, 3L))
  expect_equal(linop$data, x@id)
  expect_length(canon[[2]], 0)
})

## @cvxpy NONE
test_that("Variable print", {
  x <- Variable(3, name = "z")
  output <- capture.output(print(x))
  expect_true(grepl("Variable", output))
  expect_true(grepl("z", output))
})

# ── Constant tests ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("Constant from scalar", {
  c1 <- Constant(5)
  expect_equal(c1@shape, c(1L, 1L))
  expect_equal(value(c1), matrix(5))
  expect_true(is_constant(c1))
})

## @cvxpy NONE
test_that("Constant from vector", {
  c1 <- Constant(1:3)
  expect_equal(c1@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("Constant from matrix", {
  m <- matrix(1:6, 2, 3)
  c1 <- Constant(m)
  expect_equal(c1@shape, c(2L, 3L))
  expect_equal(value(c1), m)
})

## @cvxpy NONE
test_that("Constant from sparse matrix", {
  sm <- Matrix::sparseMatrix(i = c(1, 2), j = c(1, 3), x = c(1.0, 2.0),
                             dims = c(3, 3))
  c1 <- Constant(sm)
  expect_equal(c1@shape, c(3L, 3L))
  expect_true(c1@.sparse)
})

## @cvxpy NONE
test_that("Constant sign: positive", {
  c1 <- Constant(matrix(c(1, 2, 3), 1, 3))
  expect_true(is_nonneg(c1))
  expect_false(is_nonpos(c1))
})

## @cvxpy NONE
test_that("Constant sign: negative", {
  c1 <- Constant(matrix(c(-1, -2, -3), 1, 3))
  expect_false(is_nonneg(c1))
  expect_true(is_nonpos(c1))
})

## @cvxpy NONE
test_that("Constant sign: zero", {
  c1 <- Constant(0)
  expect_true(is_nonneg(c1))
  expect_true(is_nonpos(c1))
  expect_true(is_zero(c1))
})

## @cvxpy NONE
test_that("Constant sign: mixed", {
  c1 <- Constant(c(-1, 1))
  expect_false(is_nonneg(c1))
  expect_false(is_nonpos(c1))
  expect_false(is_zero(c1))
})

## @cvxpy NONE
test_that("Constant curvature is always CONSTANT", {
  c1 <- Constant(5)
  expect_true(is_constant(c1))
  expect_true(is_affine(c1))
  expect_true(is_convex(c1))
  expect_true(is_concave(c1))
  expect_true(is_dcp(c1))
  expect_equal(expr_curvature(c1), "CONSTANT")
})

## @cvxpy NONE
test_that("Constant is_symmetric: scalar", {
  expect_true(is_symmetric(Constant(5)))
})

## @cvxpy NONE
test_that("Constant is_symmetric: symmetric matrix", {
  m <- matrix(c(1, 2, 2, 3), 2, 2)
  expect_true(is_symmetric(Constant(m)))
})

## @cvxpy NONE
test_that("Constant is_symmetric: non-symmetric", {
  m <- matrix(c(1, 2, 3, 4), 2, 2)
  expect_false(is_symmetric(Constant(m)))
})

## @cvxpy NONE
test_that("Constant is_symmetric: non-square", {
  m <- matrix(1:6, 2, 3)
  expect_false(is_symmetric(Constant(m)))
})

## @cvxpy NONE
test_that("Constant PSD: positive definite matrix", {
  m <- matrix(c(2, 1, 1, 2), 2, 2)
  expect_true(is_psd(Constant(m)))
})

## @cvxpy NONE
test_that("Constant NSD: negative definite matrix", {
  m <- matrix(c(-2, -1, -1, -2), 2, 2)
  expect_true(is_nsd(Constant(m)))
})

## @cvxpy NONE
test_that("Constant PSD: indefinite matrix", {
  m <- matrix(c(1, 0, 0, -1), 2, 2)
  expect_false(is_psd(Constant(m)))
  expect_false(is_nsd(Constant(m)))
})

## @cvxpy NONE
test_that("Constant PSD: scalar positive", {
  expect_true(is_psd(Constant(5)))
  expect_false(is_nsd(Constant(5)))
})

## @cvxpy NONE
test_that("Constant PSD: scalar negative", {
  expect_false(is_psd(Constant(-5)))
  expect_true(is_nsd(Constant(-5)))
})

## @cvxpy NONE
test_that("Constant immutable (value<- errors)", {
  c1 <- Constant(5)
  expect_error(value(c1) <- 10, "Cannot set the value of a.*Constant")
})

## @cvxpy NONE
test_that("Constant constants() returns self", {
  c1 <- Constant(5)
  consts <- constants(c1)
  expect_length(consts, 1)
  expect_identical(consts[[1]]@id, c1@id)
})

## @cvxpy NONE
test_that("Constant variables() returns empty", {
  c1 <- Constant(5)
  expect_length(variables(c1), 0)
})

## @cvxpy NONE
test_that("Constant grad returns empty", {
  c1 <- Constant(5)
  expect_length(grad(c1), 0)
})

## @cvxpy NONE
test_that("Constant canonicalize: scalar", {
  c1 <- Constant(42)
  canon <- canonicalize(c1)
  linop <- canon[[1]]
  expect_equal(linop$type, "scalar_const")
  expect_equal(linop$data, 42)
  expect_equal(linop$shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("Constant canonicalize: matrix", {
  m <- matrix(1:6, 2, 3)
  c1 <- Constant(m)
  canon <- canonicalize(c1)
  linop <- canon[[1]]
  expect_equal(linop$type, "dense_const")
  expect_equal(linop$shape, c(2L, 3L))
})

## @cvxpy NONE
test_that("Constant canonicalize: sparse", {
  sm <- Matrix::sparseMatrix(i = c(1, 2), j = c(1, 3), x = c(1.0, 2.0),
                             dims = c(3, 3))
  c1 <- Constant(sm)
  canon <- canonicalize(c1)
  linop <- canon[[1]]
  expect_equal(linop$type, "sparse_const")
})

## @cvxpy NONE
test_that("Constant print", {
  c1 <- Constant(5)
  output <- capture.output(print(c1))
  expect_true(grepl("Constant", output))
})

## @cvxpy NONE
test_that("Constant expr_name: scalar", {
  c1 <- Constant(42)
  expect_equal(expr_name(c1), "42")
})

## @cvxpy NONE
test_that("Constant expr_name: matrix", {
  m <- matrix(1:6, 2, 3)
  c1 <- Constant(m)
  expect_true(grepl("2x3", expr_name(c1)))
})

## @cvxpy NONE
test_that("Constant expr_name: custom", {
  c1 <- Constant(5, name = "five")
  expect_equal(expr_name(c1), "five")
})

## @cvxpy NONE
test_that("Constant from logical converts to numeric", {
  c1 <- Constant(TRUE)
  expect_true(is.numeric(value(c1)))
  expect_equal(value(c1), matrix(1))
})

# ── Parameter tests ───────────────────────────────────────────────────

## @cvxpy NONE
test_that("Parameter construction", {
  CVXR:::reset_expr_id()
  p <- Parameter(c(2, 3))
  expect_equal(p@shape, c(2L, 3L))
  expect_equal(expr_name(p), "param1")
  expect_true(is_constant(p))
})

## @cvxpy NONE
test_that("Parameter scalar shape normalization", {
  p <- Parameter(5)
  expect_equal(p@shape, c(5L, 1L))
})

## @cvxpy NONE
test_that("Parameter custom name", {
  p <- Parameter(c(2, 2), name = "alpha")
  expect_equal(expr_name(p), "alpha")
})

## @cvxpy NONE
test_that("Parameter value assignment", {
  p <- Parameter(c(2, 3))
  expect_null(value(p))
  value(p) <- matrix(1:6, 2, 3)
  expect_equal(sum(value(p)), 21)
})

## @cvxpy NONE
test_that("Parameter value: wrong shape errors", {
  p <- Parameter(c(2, 3))
  expect_error(value(p) <- matrix(1:4, 2, 2), "Invalid dimensions")
})

## @cvxpy NONE
test_that("Parameter nonneg: valid assignment", {
  p <- Parameter(c(2, 2), nonneg = TRUE)
  expect_true(is_nonneg(p))
  value(p) <- matrix(c(1, 2, 3, 4), 2, 2)
  expect_equal(sum(value(p)), 10)
})

## @cvxpy NONE
test_that("Parameter nonneg: invalid assignment errors", {
  p <- Parameter(c(2, 2), nonneg = TRUE)
  expect_error(value(p) <- matrix(c(-1, 2, 3, 4), 2, 2), "nonnegative")
})

## @cvxpy NONE
test_that("Parameter nonpos: valid assignment", {
  p <- Parameter(c(2, 2), nonpos = TRUE)
  expect_true(is_nonpos(p))
  value(p) <- matrix(c(-1, -2, -3, -4), 2, 2)
  expect_equal(sum(value(p)), -10)
})

## @cvxpy NONE
test_that("Parameter nonpos: invalid assignment errors", {
  p <- Parameter(c(2, 2), nonpos = TRUE)
  expect_error(value(p) <- matrix(c(1, -2, -3, -4), 2, 2), "nonpositive")
})

## @cvxpy NONE
test_that("Parameter cache invalidation on value<-", {
  p <- Parameter(c(2, 2))
  ## Access a cached method to populate cache
  is_constant(p)
  ## Now set value — should clear cache
  value(p) <- matrix(1:4, 2, 2)
  ## Cache should be cleared
  expect_false(CVXR:::cache_has(p, "is_constant"))
})

## @cvxpy NONE
test_that("Parameter parameters() returns self", {
  p <- Parameter(c(2, 2))
  params <- parameters(p)
  expect_length(params, 1)
  expect_identical(params[[1]]@id, p@id)
})

## @cvxpy NONE
test_that("Parameter variables() returns empty", {
  p <- Parameter(c(2, 2))
  expect_length(variables(p), 0)
})

## @cvxpy NONE
test_that("Parameter constants() returns empty", {
  p <- Parameter(c(2, 2))
  expect_length(constants(p), 0)
})

## @cvxpy NONE
test_that("Parameter grad returns empty", {
  p <- Parameter(c(2, 2))
  expect_length(grad(p), 0)
})

## @cvxpy NONE
test_that("Parameter canonicalize", {
  p <- Parameter(c(2, 3))
  canon <- canonicalize(p)
  linop <- canon[[1]]
  expect_equal(linop$type, "param")
  expect_equal(linop$shape, c(2L, 3L))
  expect_equal(linop$data, p@id)
})

## @cvxpy NONE
test_that("Parameter PSD attribute", {
  p <- Parameter(c(3, 3), PSD = TRUE)
  expect_true(is_psd(p))
  expect_true(is_symmetric(p))
})

## @cvxpy NONE
test_that("Parameter get_data", {
  p <- Parameter(c(2, 3), name = "beta")
  data <- get_data(p)
  expect_length(data, 5)
  expect_equal(data[[1]], c(2L, 3L))
  expect_equal(data[[2]], "beta")
})

## @cvxpy NONE
test_that("Parameter print", {
  p <- Parameter(c(2, 3), name = "gamma")
  output <- capture.output(print(p))
  expect_true(grepl("Parameter", output))
  expect_true(grepl("gamma", output))
})

## @cvxpy NONE
test_that("Parameter value<- with NULL clears", {
  p <- Parameter(c(2, 2))
  value(p) <- matrix(1:4, 2, 2)
  expect_false(is.null(value(p)))
  value(p) <- NULL
  expect_null(value(p))
})

# ── Expression base class tests ───────────────────────────────────────

## @cvxpy NONE
test_that("as_expr promotes scalar to Constant", {
  e <- as_expr(5)
  expect_s3_class(e, "CVXR::Constant")
  expect_equal(value(e), matrix(5))
})

## @cvxpy NONE
test_that("as_expr promotes matrix to Constant", {
  m <- matrix(1:4, 2, 2)
  e <- as_expr(m)
  expect_s3_class(e, "CVXR::Constant")
  expect_equal(value(e), m)
})

## @cvxpy NONE
test_that("as_expr passes through Expression", {
  x <- Variable(3)
  e <- as_expr(x)
  expect_identical(e, x)
})

## @cvxpy NONE
test_that("as_expr errors on invalid type", {
  expect_error(as_expr("hello"), "Cannot convert")
})

## @cvxpy NONE
test_that("expr_curvature strings", {
  x <- Variable(3)
  expect_equal(expr_curvature(x), "AFFINE")
  c1 <- Constant(5)
  expect_equal(expr_curvature(c1), "CONSTANT")
})

## @cvxpy NONE
test_that("caching works (repeated calls return same result)", {
  x <- Variable(3)
  r1 <- is_affine(x)
  r2 <- is_affine(x)
  expect_identical(r1, r2)
  expect_true(CVXR:::cache_has(x, "is_affine"))
})

# ── Class hierarchy tests ────────────────────────────────────────────

## @cvxpy NONE
test_that("Variable inherits from Leaf, Expression, Canonical", {
  x <- Variable(3)
  expect_s3_class(x, "CVXR::Variable")
  expect_s3_class(x, "CVXR::Leaf")
  expect_s3_class(x, "CVXR::Expression")
  expect_s3_class(x, "CVXR::Canonical")
})

## @cvxpy NONE
test_that("Constant inherits from Leaf, Expression, Canonical", {
  c1 <- Constant(5)
  expect_s3_class(c1, "CVXR::Constant")
  expect_s3_class(c1, "CVXR::Leaf")
  expect_s3_class(c1, "CVXR::Expression")
  expect_s3_class(c1, "CVXR::Canonical")
})

## @cvxpy NONE
test_that("Parameter inherits from Leaf, Expression, Canonical", {
  p <- Parameter(c(2, 2))
  expect_s3_class(p, "CVXR::Parameter")
  expect_s3_class(p, "CVXR::Leaf")
  expect_s3_class(p, "CVXR::Expression")
  expect_s3_class(p, "CVXR::Canonical")
})

# ── Shape utility tests ──────────────────────────────────────────────

## @cvxpy NONE
test_that("expr_size computes product of shape dims", {
  x <- Variable(c(3, 4))
  expect_equal(expr_size(x), 12L)
})

## @cvxpy NONE
test_that("expr_ndim returns number of dimensions", {
  x <- Variable(c(3, 4))
  expect_equal(expr_ndim(x), 2L)
})

## @cvxpy NONE
test_that("expr_is_scalar", {
  expect_true(expr_is_scalar(Variable(c(1, 1))))
  expect_false(expr_is_scalar(Variable(3)))
})

## @cvxpy NONE
test_that("expr_is_vector", {
  expect_true(expr_is_vector(Variable(3)))
  expect_true(expr_is_vector(Variable(c(1, 5))))
  expect_false(expr_is_vector(Variable(c(2, 3))))
})

## @cvxpy NONE
test_that("expr_is_matrix", {
  expect_true(expr_is_matrix(Variable(c(2, 3))))
  expect_false(expr_is_matrix(Variable(3)))
  expect_false(expr_is_matrix(Variable(c(1, 3))))
})

# ── validate_shape tests ─────────────────────────────────────────────

## @cvxpy NONE
test_that("validate_shape: NULL becomes c(1, 1)", {
  expect_equal(CVXR:::validate_shape(NULL), c(1L, 1L))
})

## @cvxpy NONE
test_that("validate_shape: scalar becomes c(n, 1)", {
  expect_equal(CVXR:::validate_shape(3), c(3L, 1L))
})

## @cvxpy NONE
test_that("validate_shape: 2D passes through", {
  expect_equal(CVXR:::validate_shape(c(2, 3)), c(2L, 3L))
})

## @cvxpy NONE
test_that("validate_shape: rejects 3D", {
  expect_error(CVXR:::validate_shape(c(2, 3, 4)), "length 1 or 2")
})

## @cvxpy NONE
test_that("validate_shape: rejects negative", {
  expect_error(CVXR:::validate_shape(c(-1, 2)), "positive")
})

## @cvxpy NONE
test_that("validate_shape: rejects zero dimensions", {
  expect_error(CVXR:::validate_shape(c(0, 2)), "positive")
  expect_error(CVXR:::validate_shape(c(3, 0)), "positive")
})

# ── Leaf project tests ───────────────────────────────────────────────

## @cvxpy NONE
test_that("project: nonneg clips negatives", {
  x <- Variable(3, nonneg = TRUE)
  val <- matrix(c(-1, 2, -3))
  result <- project(x, val)
  expect_equal(as.numeric(result), c(0, 2, 0))
})

## @cvxpy NONE
test_that("project: nonpos clips positives", {
  x <- Variable(3, nonpos = TRUE)
  val <- matrix(c(-1, 2, -3))
  result <- project(x, val)
  expect_equal(as.numeric(result), c(-1, 0, -3))
})

## @cvxpy NONE
test_that("project: symmetric averages", {
  x <- Variable(c(2, 2), symmetric = TRUE)
  val <- matrix(c(1, 2, 3, 4), 2, 2)
  result <- project(x, val)
  expect_equal(result, (val + t(val)) / 2)
})

## @cvxpy NONE
test_that("project: PSD clamps negative eigenvalues", {
  x <- Variable(c(2, 2), PSD = TRUE)
  val <- matrix(c(1, 0, 0, -1), 2, 2)
  result <- project(x, val)
  ev <- eigen(result, symmetric = TRUE)$values
  expect_true(all(ev >= -1e-10))
})

## @cvxpy NONE
test_that("project: boolean rounds", {
  x <- Variable(3, boolean = TRUE)
  val <- matrix(c(0.3, 0.7, -0.2))
  result <- project(x, val)
  expect_equal(as.numeric(result), c(0, 1, 0))
})

## @cvxpy NONE
test_that("project: integer rounds", {
  x <- Variable(3, integer = TRUE)
  val <- matrix(c(1.3, 2.7, -0.8))
  result <- project(x, val)
  expect_equal(as.numeric(result), c(1, 3, -1))
})

# ── Interface utility tests ──────────────────────────────────────────

## @cvxpy NONE
test_that("intf_convert: scalar to 1x1 matrix", {
  result <- CVXR:::intf_convert(42)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(1L, 1L))
  expect_equal(result[1, 1], 42)
})

## @cvxpy NONE
test_that("intf_convert: vector to column matrix", {
  result <- CVXR:::intf_convert(1:3)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3L, 1L))
})

## @cvxpy NONE
test_that("intf_convert: matrix passes through", {
  m <- matrix(1:6, 2, 3)
  result <- CVXR:::intf_convert(m)
  expect_equal(result, m)
})

## @cvxpy NONE
test_that("intf_convert: logical to numeric", {
  result <- CVXR:::intf_convert(c(TRUE, FALSE))
  expect_true(is.numeric(result))
  expect_equal(dim(result), c(2L, 1L))
})

## @cvxpy NONE
test_that("intf_shape: matrix", {
  expect_equal(CVXR:::intf_shape(matrix(0, 3, 4)), c(3L, 4L))
})

## @cvxpy NONE
test_that("intf_shape: scalar", {
  expect_equal(CVXR:::intf_shape(5), c(1L, 1L))
})

## @cvxpy NONE
test_that("intf_sign: positive", {
  sgn <- CVXR:::intf_sign(matrix(c(1, 2, 3)))
  expect_true(sgn$is_nonneg)
  expect_false(sgn$is_nonpos)
})

## @cvxpy NONE
test_that("intf_sign: mixed", {
  sgn <- CVXR:::intf_sign(matrix(c(-1, 2)))
  expect_false(sgn$is_nonneg)
  expect_false(sgn$is_nonpos)
})

## @cvxpy NONE
test_that("intf_is_sparse: detects sparse", {
  sm <- Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(2, 2))
  expect_true(CVXR:::intf_is_sparse(sm))
  expect_false(CVXR:::intf_is_sparse(matrix(0, 2, 2)))
})

## @cvxpy NONE
test_that("intf_is_hermitian: symmetric", {
  m <- matrix(c(1, 2, 2, 3), 2, 2)
  result <- CVXR:::intf_is_hermitian(m)
  expect_true(result$is_symmetric)
})

## @cvxpy NONE
test_that("intf_is_hermitian: non-symmetric", {
  m <- matrix(c(1, 2, 3, 4), 2, 2)
  result <- CVXR:::intf_is_hermitian(m)
  expect_false(result$is_symmetric)
})

# ── LinOp utility tests ──────────────────────────────────────────────

## @cvxpy NONE
test_that("create_var produces correct LinOp", {
  lo <- CVXR:::create_var(c(3L, 1L), 42L)
  expect_equal(lo$type, "variable")
  expect_equal(lo$shape, c(3L, 1L))
  expect_equal(lo$data, 42L)
  expect_equal(class(lo), "LinOp_R")
})

## @cvxpy NONE
test_that("create_param produces correct LinOp", {
  lo <- CVXR:::create_param(c(2L, 3L), 99L)
  expect_equal(lo$type, "param")
  expect_equal(lo$shape, c(2L, 3L))
  expect_equal(lo$data, 99L)
})

## @cvxpy NONE
test_that("create_const: scalar", {
  lo <- CVXR:::create_const(5, c(1L, 1L), FALSE)
  expect_equal(lo$type, "scalar_const")
  expect_equal(lo$data, 5)
})

## @cvxpy NONE
test_that("create_const: dense matrix", {
  m <- matrix(1:6, 2, 3)
  lo <- CVXR:::create_const(m, c(2L, 3L), FALSE)
  expect_equal(lo$type, "dense_const")
  expect_equal(lo$shape, c(2L, 3L))
})

## @cvxpy NONE
test_that("create_const: sparse matrix", {
  sm <- Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(3, 3))
  lo <- CVXR:::create_const(sm, c(3L, 3L), TRUE)
  expect_equal(lo$type, "sparse_const")
})

# ── Phase 1 Review Fix Tests ─────────────────────────────────────────

# B1: Sign string constants match CVXPY
## @cvxpy NONE
test_that("sign constants match CVXPY (NONNEGATIVE/NONPOSITIVE)", {
  expect_equal(NONNEG_SIGN, "NONNEGATIVE")
  expect_equal(NONPOS_SIGN, "NONPOSITIVE")
  expect_equal(ZERO_SIGN, "ZERO")
  expect_equal(UNKNOWN_SIGN, "UNKNOWN")
})

## @cvxpy NONE
test_that("expr_sign_str uses correct CVXPY sign strings", {
  c_pos <- Constant(5)
  expect_equal(CVXR:::expr_sign_str(c_pos), "NONNEGATIVE")
  c_neg <- Constant(-3)
  expect_equal(CVXR:::expr_sign_str(c_neg), "NONPOSITIVE")
  c_zero <- Constant(0)
  expect_equal(CVXR:::expr_sign_str(c_zero), "ZERO")
  c_mixed <- Constant(matrix(c(-1, 1)))
  expect_equal(CVXR:::expr_sign_str(c_mixed), "UNKNOWN")
})

# B2: canonical_form lazy caching
## @cvxpy NONE
test_that("canonical_form caches canonicalize result", {
  x <- Variable(2)
  cf1 <- canonical_form(x)
  cf2 <- canonical_form(x)
  expect_identical(cf1, cf2)
  expect_equal(cf1[[1]]$type, "variable")
  expect_equal(cf1[[2]], list())
})

## @cvxpy NONE
test_that("canonical_form works for Constant", {
  c <- Constant(42)
  cf <- canonical_form(c)
  expect_equal(cf[[1]]$type, "scalar_const")
  expect_equal(cf[[1]]$data, 42)
})

## @cvxpy NONE
test_that("canonical_form works for Parameter", {
  p <- Parameter(c(2, 3))
  cf <- canonical_form(p)
  expect_equal(cf[[1]]$type, "param")
})

# B3: copy/tree_copy stubs exist
## @cvxpy NONE
test_that("expr_copy returns self for leaf nodes", {
  x <- Variable(3)
  copied <- expr_copy(x)
  expect_identical(copied, x)
})

## @cvxpy NONE
test_that("tree_copy returns self for leaf nodes", {
  x <- Variable(3)
  copied <- tree_copy(x)
  expect_identical(copied, x)
})

# H2: intf_sign handles NaN and empty
## @cvxpy NONE
test_that("intf_sign: NaN returns unknown sign", {
  result <- CVXR:::intf_sign(matrix(c(1, NaN, 3)))
  expect_false(result$is_nonneg)
  expect_false(result$is_nonpos)
})

## @cvxpy NONE
test_that("intf_sign: empty matrix returns both TRUE", {
  m <- matrix(numeric(0), nrow = 0, ncol = 0)
  result <- CVXR:::intf_sign(m)
  expect_true(result$is_nonneg)
  expect_true(result$is_nonpos)
})

## @cvxpy NONE
test_that("intf_sign: zero sparse matrix", {
  sm <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                              dims = c(3, 3))
  result <- CVXR:::intf_sign(sm)
  expect_true(result$is_nonneg)
  expect_true(result$is_nonpos)
})

# H6: project with bounds
## @cvxpy NONE
test_that("project: bounds clamping", {
  x <- Variable(3, bounds = list(-1, 2))
  val <- matrix(c(-5, 0, 10))
  result <- project(x, val)
  expect_equal(as.numeric(result), c(-1, 0, 2))
})

# H8: Parameter value<- uses shared validation
## @cvxpy NONE
test_that("Parameter value<- validates like Leaf", {
  p <- Parameter(c(2, 2), nonneg = TRUE)
  expect_error(value(p) <- matrix(c(-1, 2, 3, 4), 2, 2), "nonnegative")
  value(p) <- matrix(c(1, 2, 3, 4), 2, 2)
  expect_equal(as.numeric(value(p)), c(1, 2, 3, 4))
})

# H9: PSD_NSD_PROJECTION_TOL matches CVXPY (1e-8)
## @cvxpy NONE
test_that("PSD_NSD_PROJECTION_TOL is 1e-8", {
  expect_equal(CVXR:::PSD_NSD_PROJECTION_TOL, 1e-8)
})

# New generics: is_hermitian, is_real, is_log_log_*
## @cvxpy NONE
test_that("is_hermitian defaults to is_symmetric", {
  x <- Variable(c(2, 2), symmetric = TRUE)
  expect_true(is_hermitian(x))
  y <- Variable(c(2, 3))
  expect_false(is_hermitian(y))
})

## @cvxpy NONE
test_that("is_real default is TRUE for non-complex", {
  x <- Variable(3)
  expect_true(is_real(x))
})

## @cvxpy NONE
test_that("is_log_log_* defaults to FALSE", {
  x <- Variable(3)
  expect_false(is_log_log_convex(x))
  expect_false(is_log_log_concave(x))
  expect_false(is_log_log_affine(x))
})

# M6: unique_list uses O(1) lookup
## @cvxpy NONE
test_that("unique_list deduplicates by id efficiently", {
  x1 <- Variable(2)
  x2 <- Variable(3)
  lst <- list(x1, x2, x1, x2, x1)
  result <- CVXR:::unique_list(lst)
  expect_length(result, 2)
  expect_equal(result[[1]]@id, x1@id)
  expect_equal(result[[2]]@id, x2@id)
})
