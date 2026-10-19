## Masks of base/stats/Matrix names on PLAIN R objects: CVXR is not
## involved, so the result must be identical() to the original function
## called exactly as written. CVXR-only arguments on plain objects error.

library(testthat)
library(CVXR)

## ═══════════════════════════════════════════════════════════════════
## diag
## ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("diag mask on plain objects is identical to Matrix::diag as written", {
  A <- matrix(1:4, 2, dimnames = list(c("a", "b"), c("a", "b")))
  v <- c(1, 2, 3)
  S <- Matrix::Matrix(c(1, 0, 0, 2), 2, sparse = TRUE)
  ## base diag() with no arguments errors on the missing nrow; so must the mask
  expect_error(CVXR::diag(), "nrow")
  expect_identical(CVXR::diag(3), Matrix::diag(3))
  expect_identical(CVXR::diag(nrow = 3), Matrix::diag(nrow = 3))
  expect_identical(CVXR::diag(2, 3, 4), Matrix::diag(2, 3, 4))
  expect_identical(CVXR::diag(v), Matrix::diag(v))
  expect_identical(CVXR::diag(v, nrow = 4), Matrix::diag(v, nrow = 4))
  expect_identical(CVXR::diag(A), Matrix::diag(A))
  expect_identical(CVXR::diag(A, names = FALSE), Matrix::diag(A, names = FALSE))
  expect_identical(CVXR::diag(S), Matrix::diag(S))
})

## @cvxpy NONE
test_that("diag mask: k on a plain object is an error, not ignored", {
  A <- matrix(1:4, 2)
  expect_error(CVXR::diag(A, k = 1), "k")
  expect_error(CVXR::diag(c(1, 2), k = 0), "k")
})

## @cvxpy NONE
test_that("diag mask on expressions keeps the CVXPY k default of 0", {
  X <- Variable(c(3, 3))
  x <- Variable(3)
  expect_equal(dim(CVXR::diag(X)), c(3L, 1L))
  expect_equal(dim(CVXR::diag(X, k = 1)), c(2L, 1L))
  expect_equal(dim(CVXR::diag(x)), c(3L, 3L))
  expect_equal(dim(CVXR::diag(x, k = 1)), c(4L, 4L))
})

## diag.py:27 diag(expr, k=0) has no nrow/ncol/names; they were ignored.
## @cvxpy NONE
test_that("diag mask: nrow, ncol and names on an expression are errors", {
  X <- Variable(c(3, 3))
  x <- Variable(3)
  expect_error(CVXR::diag(x, 2), "nrow")
  expect_error(CVXR::diag(x, nrow = 4), "nrow")
  expect_error(CVXR::diag(x, ncol = 4), "ncol")
  expect_error(CVXR::diag(X, names = FALSE), "names")
  expect_error(CVXR::diag(x, 4, 4, k = 1), "nrow")
})

## ═══════════════════════════════════════════════════════════════════
## norm
## ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("norm mask on plain objects is identical to Matrix::norm as written", {
  A <- matrix(c(1, -2, 3, 4), 2)
  S <- Matrix::Matrix(c(1, 0, -3, 2), 2, sparse = TRUE)
  ## no type: the base/Matrix default ("O") applies, not CVXR's "2"
  expect_identical(CVXR::norm(A), Matrix::norm(A))
  expect_identical(CVXR::norm(S), Matrix::norm(S))
  for (tp in c("O", "I", "F", "M", "2", "e", "1"))
    expect_identical(CVXR::norm(A, tp), Matrix::norm(A, tp))
  expect_identical(CVXR::norm(A, type = "F"), Matrix::norm(A, type = "F"))
  expect_identical(CVXR::norm(S, "I"), Matrix::norm(S, "I"))
})

## @cvxpy NONE
test_that("norm mask on expressions keeps the spectral default", {
  X <- Variable(c(2, 2))
  expect_identical(class(CVXR::norm(X)), class(cvxr_norm(X, p = 2)))
  expect_identical(class(CVXR::norm(X, "O")), class(cvxr_norm(X, p = 1)))
})

## ═══════════════════════════════════════════════════════════════════
## convolve
## ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("convolve mask on plain objects is identical to stats::convolve as written", {
  v <- c(1, 2, 3, 4)
  w <- c(1, -1, 2, 0.5)
  k <- c(1, 1)
  ## no type: stats' circular default, not numpy's open convolution
  expect_identical(CVXR::convolve(v, w), stats::convolve(v, w))
  expect_identical(CVXR::convolve(v, w, conj = FALSE),
                   stats::convolve(v, w, conj = FALSE))
  expect_identical(CVXR::convolve(v, k, type = "open"),
                   stats::convolve(v, k, type = "open"))
  expect_identical(CVXR::convolve(v, k, type = "filter"),
                   stats::convolve(v, k, type = "filter"))
  expect_identical(CVXR::convolve(y = k, x = v, type = "o"),
                   stats::convolve(y = k, x = v, type = "o"))
})

## @cvxpy NONE
test_that("convolve mask on expressions builds the atom; conj/type are stats-only", {
  x <- Variable(3)
  expect_s3_class(CVXR::convolve(c(1, 1), x), "CVXR::Convolve")
  expect_equal(dim(CVXR::convolve(c(1, 1), x)), c(4L, 1L))
  expect_error(CVXR::convolve(c(1, 1), x, type = "open"), "only for plain")
  expect_error(CVXR::convolve(c(1, 1), x, conj = FALSE), "only for plain")
})

## ═══════════════════════════════════════════════════════════════════
## power: NOT a mask. stats::power() is the unrelated glm() link
## constructor; power() and cvxr_power() are CVXPY's power(x, p) for any
## argument, as in CVXR 1.9.2.
## ═══════════════════════════════════════════════════════════════════

## power.py:31 power(x, p, max_denom=1024, approx=True). The scpi package
## builds constraint bounds with CVXR::power(QQ, 2) on a plain number.
## @cvxpy NONE
test_that("power() and cvxr_power() build the atom for plain numbers too", {
  expect_identical(CVXR::power, CVXR::cvxr_power)
  pq <- CVXR::power(3, 2)
  expect_s3_class(pq, "CVXR::Expression")
  expect_equal(as.numeric(value(pq)), 9)
  expect_equal(as.numeric(value(cvxr_power(2, p = 3))), 8)
  x <- Variable(3)
  expect_s3_class(CVXR::power(x, 2), "CVXR::PowerApprox")
  expect_s3_class(CVXR::power(x, p = 2), "CVXR::PowerApprox")
  expect_s3_class(CVXR::power(x, 2, approx = FALSE), "CVXR::Power")
  expect_s3_class(CVXR::power(2, x), "CVXR::Expression")
  expect_s3_class(square(3), "CVXR::Expression")
  expect_s3_class(inv_pos(2), "CVXR::Expression")
  ## the scpi form: a bound built from a plain number
  prob <- Problem(Minimize(sum(x)), list(sum_squares(x) <= CVXR::power(3, 2)))
  expect_equal(psolve(prob, solver = "CLARABEL"), -3 * sqrt(3), tolerance = 1e-6)
})

## ═══════════════════════════════════════════════════════════════════
## sd, var
## ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("sd/var masks on plain objects are identical to stats as written", {
  v <- c(2, 3, 4.5, 5, 1, NA)
  m <- matrix(c(2, 3, 4.5, 5, 1, 7), 3, 2)
  expect_identical(sd(v, na.rm = TRUE), stats::sd(v, na.rm = TRUE))
  expect_identical(sd(v, TRUE), stats::sd(v, TRUE))
  expect_identical(var(m), stats::var(m))
  expect_identical(var(m[, 1], m[, 2]), stats::var(m[, 1], m[, 2]))
  expect_identical(var(v, use = "complete.obs"), stats::var(v, use = "complete.obs"))
})

## stats.py:42/57 std/var(x, axis=None, keepdims=False, ddof=0); any other
## keyword is a TypeError there.
## @cvxpy NONE
test_that("sd/var masks on expressions are cvxr_std/cvxr_var (ddof 0)", {
  m <- matrix(c(2, 3, 4.5, 5, 1, 7), 2, 3)
  n <- length(m)
  expect_equal(as.numeric(value(sd(Constant(m)))), sqrt(sum((m - mean(m))^2) / n))
  expect_equal(as.numeric(value(var(Constant(m)))), sum((m - mean(m))^2) / n)
  expect_equal(as.numeric(value(var(Constant(m), ddof = 1))), stats::var(as.vector(m)))
  x <- Variable(c(2, 3))
  expect_error(sd(x, na.rm = TRUE), "cvxr_std")
  expect_error(var(x, na.rm = TRUE), "cvxr_var")
  expect_error(var(x, use = "all.obs"), "cvxr_var")
  expect_error(var(x, y = x), "cvxr_var")
  ## a positional second argument is cvxr_var's axis
  expect_error(var(x, x), "axis")
  expect_error(var(x, m), "axis")
  expect_error(sd(x, c(1, 2)), "axis")
})

## ═══════════════════════════════════════════════════════════════════
## outer
## ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("outer mask on plain objects is identical to base::outer as written", {
  a <- 1:3
  b <- c(2, 5)
  expect_identical(outer(a, b), base::outer(a, b))
  expect_identical(outer(a, b, "+"), base::outer(a, b, "+"))
  expect_identical(outer(a, b, FUN = paste, sep = "-"),
                   base::outer(a, b, FUN = paste, sep = "-"))
})

## binary_operators.py:682 outer(x, y) has no FUN; on an expression only the
## product is meaningful, so any other FUN or extra argument is an error.
## @cvxpy NONE
test_that("outer mask on expressions is the product; another FUN is an error", {
  a <- c(1, 2, 3)
  b <- c(2, 5)
  expect_equal(value(outer(Constant(a), b)), base::outer(a, b))
  expect_equal(value(outer(Constant(a), b, "*")), base::outer(a, b))
  expect_equal(value(outer(a, Constant(b), FUN = `*`)), base::outer(a, b))
  x <- Variable(3)
  expect_error(outer(x, b, "+"), "FUN")
  expect_error(outer(x, b, FUN = "-"), "FUN")
  expect_error(outer(b, x, function(u, v) u * v), "FUN")
  expect_error(outer(x, b, "*", 2), "FUN")
  expect_error(outer(x, b, sep = "-"), "FUN")
})
