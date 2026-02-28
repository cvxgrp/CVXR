## Tests for Phase 7a: Wrappers, Utilities, Quick Wins

# ═══════════════════════════════════════════════════════════════════
# diff() tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("cvxr_diff works on column vector", {
  x <- Variable(c(5L, 1L))
  d <- cvxr_diff(x)
  ## Shape: 5 rows - 1 = 4 rows

  expect_equal(d@shape, c(4L, 1L))
})

## @cvxpy NONE
test_that("cvxr_diff works on matrix, axis=0 (rows)", {
  x <- Variable(c(4L, 3L))
  d <- cvxr_diff(x, axis = 2L)
  expect_equal(d@shape, c(3L, 3L))
})

## @cvxpy NONE
test_that("cvxr_diff works on matrix, axis=1 (columns)", {
  x <- Variable(c(4L, 3L))
  d <- cvxr_diff(x, axis = 1L)
  expect_equal(d@shape, c(4L, 2L))
})

## @cvxpy NONE
test_that("cvxr_diff with k=2", {
  x <- Variable(c(5L, 1L))
  d <- cvxr_diff(x, k = 2L)
  expect_equal(d@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("cvxr_diff with k=0 returns identity", {
  x <- Variable(c(5L, 1L))
  d <- cvxr_diff(x, k = 0L)
  ## Should return the same expression (no differencing)
  expect_true(S7_inherits(d, Expression))
  expect_equal(d@shape, c(5L, 1L))
})

## @cvxpy NONE
test_that("cvxr_diff validates arguments", {
  x <- Variable(c(5L, 1L))
  expect_error(cvxr_diff(x, k = -1L))
  expect_error(cvxr_diff(x, k = 5L))
  expect_error(cvxr_diff(x, axis = 3L))
})

## @cvxpy NONE
test_that("diff() S3 dispatch works on Expression", {
  x <- Variable(c(5L, 1L))
  d <- diff(x)
  expect_equal(d@shape, c(4L, 1L))
})

## @cvxpy NONE
test_that("diff() S3 dispatch with differences argument", {
  x <- Variable(c(5L, 1L))
  d <- diff(x, differences = 2L)
  expect_equal(d@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("diff() numeric evaluation", {
  x <- Variable(c(5L, 1L))
  value(x) <- matrix(c(1, 3, 6, 10, 15), 5, 1)
  d <- cvxr_diff(x)
  expect_true(S7_inherits(d, Expression))
  dval <- value(d)
  expect_equal(as.numeric(dval), c(2, 3, 4, 5), tolerance = 1e-10)
})

## @cvxpy NONE
test_that("diff() numeric evaluation axis=1", {
  x <- Variable(c(2L, 4L))
  value(x) <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), 2, 4)
  d <- cvxr_diff(x, axis = 1L)
  expect_equal(d@shape, c(2L, 3L))
  dval <- value(d)
  ## Each row: differences along columns
  ## Row 1: 3-1=2, 5-3=2, 7-5=2
  ## Row 2: 4-2=2, 6-4=2, 8-6=2
  expect_equal(as.numeric(dval), c(2, 2, 2, 2, 2, 2), tolerance = 1e-10)
})

## @cvxpy NONE
test_that("diff() is affine (DCP)", {
  x <- Variable(c(5L, 1L))
  d <- cvxr_diff(x)
  expect_true(is_affine(d))
})

# ═══════════════════════════════════════════════════════════════════
# sum_squares() tests (axis/keepdims update)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("sum_squares with axis=NULL (default)", {
  x <- Variable(c(3L, 1L))
  ss <- sum_squares(x)
  expect_equal(ss@shape, c(1L, 1L))
  expect_true(is_convex(ss))
})

## @cvxpy NONE
test_that("sum_squares with axis=2", {
  x <- Variable(c(3L, 2L))
  ss <- sum_squares(x, axis = 2L)
  ## axis=2 reduces rows (column-wise) on (3,2) → (1,2)
  expect_equal(ss@shape, c(1L, 2L))
})

## @cvxpy NONE
test_that("sum_squares DCP properties", {
  x <- Variable(c(3L, 1L))
  ss <- sum_squares(x)
  expect_true(is_convex(ss))
  expect_false(is_concave(ss))
  expect_true(is_nonneg(ss))
})

# ═══════════════════════════════════════════════════════════════════
# bmat() tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("bmat constructs block matrix from constants", {
  A <- Constant(matrix(1, 2, 2))
  B <- Constant(matrix(2, 2, 3))
  C <- Constant(matrix(3, 3, 2))
  D <- Constant(matrix(4, 3, 3))

  m <- bmat(list(list(A, B), list(C, D)))
  expect_equal(m@shape, c(5L, 5L))
})

## @cvxpy NONE
test_that("bmat with variables", {
  x <- Variable(c(2L, 2L))
  I2 <- Constant(diag(2))
  m <- bmat(list(list(x, I2), list(I2, x)))
  expect_equal(m@shape, c(4L, 4L))
  expect_true(is_affine(m))
})

## @cvxpy NONE
test_that("bmat single row", {
  A <- Variable(c(2L, 2L))
  B <- Variable(c(2L, 3L))
  m <- bmat(list(list(A, B)))
  expect_equal(m@shape, c(2L, 5L))
})

## @cvxpy NONE
test_that("bmat validates input", {
  expect_error(bmat(list()))
  expect_error(bmat(list(list())))
})

# ═══════════════════════════════════════════════════════════════════
# vec() tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("vec reshapes matrix to column vector", {
  x <- Variable(c(3L, 2L))
  v <- vec(x)
  expect_equal(v@shape, c(6L, 1L))
  expect_true(is_affine(v))
})

## @cvxpy NONE
test_that("vec of scalar", {
  x <- Variable(c(1L, 1L))
  v <- vec(x)
  expect_equal(v@shape, c(1L, 1L))
})

## @cvxpy NONE
test_that("vec numeric evaluation", {
  x <- Variable(c(2L, 3L))
  value(x) <- matrix(1:6, 2, 3)
  v <- vec(x)
  ## F-order (column-major): 1, 2, 3, 4, 5, 6
  vval <- value(v)
  expect_equal(as.numeric(vval), c(1, 2, 3, 4, 5, 6))
})

# ═══════════════════════════════════════════════════════════════════
# sum_smallest() tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("sum_smallest is negation of sum_largest of negation", {
  x <- Variable(c(5L, 1L))
  ss <- sum_smallest(x, 2)
  expect_equal(ss@shape, c(1L, 1L))
  expect_true(is_concave(ss))  ## sum_smallest is concave
})

## @cvxpy NONE
test_that("sum_smallest numeric evaluation", {
  x <- Variable(c(5L, 1L))
  value(x) <- matrix(c(10, 1, 5, 3, 7), 5, 1)
  ss <- sum_smallest(x, 2)
  ## Smallest 2 values: 1, 3 → sum = 4
  expect_equal(as.numeric(value(ss)), 4, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════
# vec_to_upper_tri() tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("vec_to_upper_tri strict=TRUE inverts upper_tri", {
  n <- 4L
  ell <- (n * (n - 1L)) %/% 2L  # 6
  z <- Variable(c(ell, 1L))
  M <- vec_to_upper_tri(z, strict = TRUE)
  expect_equal(M@shape, c(n, n))
  expect_true(is_affine(M))
})

## @cvxpy NONE
test_that("vec_to_upper_tri strict=FALSE includes diagonal", {
  n <- 3L
  ell <- (n * (n + 1L)) %/% 2L  # 6
  z <- Variable(c(ell, 1L))
  M <- vec_to_upper_tri(z, strict = FALSE)
  expect_equal(M@shape, c(n, n))
})

## @cvxpy NONE
test_that("vec_to_upper_tri numeric eval strict=TRUE", {
  n <- 3L
  ell <- (n * (n - 1L)) %/% 2L  # 3
  z <- Variable(c(ell, 1L))
  value(z) <- matrix(c(1, 2, 3), ell, 1)
  M <- vec_to_upper_tri(z, strict = TRUE)
  Mval <- value(M)
  ## Expected 3x3 strict upper-tri (row-major: (0,1)=1, (0,2)=2, (1,2)=3)
  expected <- matrix(c(0, 0, 0, 1, 0, 0, 2, 3, 0), 3, 3)
  expect_equal(as.matrix(Mval), expected, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("vec_to_upper_tri numeric eval strict=FALSE", {
  n <- 3L
  ell <- (n * (n + 1L)) %/% 2L  # 6
  z <- Variable(c(ell, 1L))
  value(z) <- matrix(c(1, 2, 3, 4, 5, 6), ell, 1)
  M <- vec_to_upper_tri(z, strict = FALSE)
  Mval <- value(M)
  ## Expected 3x3 upper-tri (row-major: (0,0)=1, (0,1)=2, (0,2)=3, (1,1)=4, (1,2)=5, (2,2)=6)
  expected <- matrix(c(1, 0, 0, 2, 4, 0, 3, 5, 6), 3, 3)
  expect_equal(as.matrix(Mval), expected, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("vec_to_upper_tri validates triangular number", {
  z <- Variable(c(5L, 1L))
  expect_error(vec_to_upper_tri(z, strict = TRUE), "not a strict triangular number")
})

# ═══════════════════════════════════════════════════════════════════
# cumsum_canon tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("cumsum_canon is registered", {
  expect_true(has_dcp_canon(Cumsum(Variable(3))))
})

## @cvxpy NONE
test_that("cumsum_canon produces correct structure for axis=0", {
  x <- Variable(c(4L, 1L))
  cs <- Cumsum(x, axis = 2L)
  result <- dcp_canonicalize(cs, cs@args)
  expect_true(S7_inherits(result[[1L]], Variable))
  expect_equal(result[[1L]]@shape, c(4L, 1L))
  expect_equal(length(result[[2L]]), 2L)  # 2 equality constraints
})

## @cvxpy NONE
test_that("cumsum_canon identity for single element", {
  x <- Variable(c(1L, 3L))
  cs <- Cumsum(x, axis = 2L)
  result <- dcp_canonicalize(cs, cs@args)
  ## Single row along axis 0 → identity (no constraints)
  expect_equal(length(result[[2L]]), 0L)
})

## @cvxpy NONE
test_that("cumsum_canon for axis=1", {
  x <- Variable(c(3L, 4L))
  cs <- Cumsum(x, axis = 1L)
  result <- dcp_canonicalize(cs, cs@args)
  expect_true(S7_inherits(result[[1L]], Variable))
  expect_equal(result[[1L]]@shape, c(3L, 4L))
  expect_equal(length(result[[2L]]), 2L)
})

## @cvxpy NONE
test_that("cumsum_canon for axis=NULL", {
  x <- Variable(c(2L, 3L))
  cs <- Cumsum(x)  # axis=NULL
  result <- dcp_canonicalize(cs, cs@args)
  ## axis=NULL: flattened to (6, 1), then cumsum axis=0
  expect_true(S7_inherits(result[[1L]], Variable))
  expect_equal(result[[1L]]@shape, c(6L, 1L))
  expect_equal(length(result[[2L]]), 2L)
})

# ═══════════════════════════════════════════════════════════════════
# cumsum solve test (end-to-end via canonicalization)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("cumsum solves end-to-end with Clarabel", {
  skip_if_not_installed("clarabel")
  ## Minimize sum(cumsum(x)) s.t. x >= 1
  ## x = (x1, x2, x3), cumsum = (x1, x1+x2, x1+x2+x3)
  ## sum(cumsum) = 3*x1 + 2*x2 + x3
  ## Minimized at x = (1, 1, 1) → 3+2+1 = 6
  x <- Variable(c(3L, 1L))
  obj <- Minimize(sum(cumsum(x)))
  prob <- Problem(obj, list(x >= 1))
  result <- solve(prob, solver = "CLARABEL")
  expect_equal(result$value, 6, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 1, 1), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("cumsum solves end-to-end with axis=0", {
  skip_if_not_installed("clarabel")
  ## x is 3x1 vector, cumsum along axis=0
  x <- Variable(c(3L, 1L))
  cs <- cumsum_axis(x, axis = 2L)
  obj <- Minimize(sum(cs))
  prob <- Problem(obj, list(x >= 1))
  result <- solve(prob, solver = "CLARABEL")
  expect_equal(result$value, 6, tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════
# dcp_canonicalize generic check
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("dcp_canonicalize is an S7 generic", {
  expect_true(inherits(dcp_canonicalize, "S7_generic"))
})

# ═══════════════════════════════════════════════════════════════════
# Integration: diff in optimization
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("diff in minimize total variation", {
  skip_if_not_installed("clarabel")
  ## Minimize ||diff(x)||_1 s.t. x[1]=0, x[5]=4
  ## This is total variation minimization → piecewise constant solution
  ## Optimal: x = (0, 1, 2, 3, 4), diff = (1, 1, 1, 1), TV = 4
  x <- Variable(c(5L, 1L))
  obj <- Minimize(norm1(cvxr_diff(x)))
  prob <- Problem(obj, list(x[1L, ] == 0, x[5L, ] == 4))
  result <- solve(prob, solver = "CLARABEL")
  expect_equal(result$value, 4, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════
# pos() / neg() already exist — verify they still work
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("pos() and neg() work correctly", {
  x <- Variable(c(3L, 1L))
  p <- pos(x)
  n <- neg(x)
  expect_true(is_convex(p))
  expect_true(is_convex(n))
})

# ═══════════════════════════════════════════════════════════════════
# square() / inv_pos() already exist — verify
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("square() and inv_pos() work correctly", {
  x <- Variable(c(2L, 1L))
  sq <- square(x)
  ip <- inv_pos(x)
  expect_true(is_convex(sq))
  expect_true(is_convex(ip))
})
