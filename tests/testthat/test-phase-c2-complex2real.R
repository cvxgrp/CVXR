## Tests for Phase C2: Complex2Real reduction
## Tests the reduction in isolation (no solver invocation)

# ── complex2real_accepts ──────────────────────────────────────────

## @cvxpy NONE
test_that("complex2real_accepts detects complex problems", {
  x <- Variable(2)
  prob_real <- Problem(Minimize(sum_entries(x)), list(x >= 0))
  expect_false(complex2real_accepts(prob_real))

  z <- Variable(2, complex = TRUE)
  prob_complex <- Problem(Minimize(sum_entries(Re(z))), list())
  expect_true(complex2real_accepts(prob_complex))
})

## @cvxpy NONE
test_that("complex2real_accepts detects complex constants", {
  x <- Variable(2)
  c_complex <- Constant(c(1 + 1i, 2 + 2i))
  prob <- Problem(Minimize(sum_entries(x)), list(x == c_complex))
  expect_true(complex2real_accepts(prob))
})

# ── Constant canon ────────────────────────────────────────────────

## @cvxpy NONE
test_that("c2r_constant_canon handles real constant", {
  c_real <- Constant(3)
  result <- c2r_constant_canon(c_real, list(), list(), NULL)
  expect_false(is.null(result[[1L]]))
  expect_null(result[[2L]])
  expect_equal(value(result[[1L]]), matrix(3))
})

## @cvxpy NONE
test_that("c2r_constant_canon handles purely imaginary constant", {
  c_imag <- Constant(3i)
  result <- c2r_constant_canon(c_imag, list(), list(), NULL)
  expect_null(result[[1L]])
  expect_false(is.null(result[[2L]]))
  expect_equal(value(result[[2L]]), matrix(3))
})

## @cvxpy NONE
test_that("c2r_constant_canon handles complex constant", {
  c_cplx <- Constant(2 + 3i)
  result <- c2r_constant_canon(c_cplx, list(), list(), NULL)
  expect_false(is.null(result[[1L]]))
  expect_false(is.null(result[[2L]]))
  expect_equal(value(result[[1L]]), matrix(2))
  expect_equal(value(result[[2L]]), matrix(3))
})

# ── Variable canon ────────────────────────────────────────────────

## @cvxpy NONE
test_that("c2r_variable_canon handles real variable", {
  x <- Variable(2)
  real2imag <- new.env(hash = TRUE, parent = emptyenv())
  result <- c2r_variable_canon(x, list(), list(), real2imag)
  expect_identical(result[[1L]], x)
  expect_null(result[[2L]])
})

## @cvxpy NONE
test_that("c2r_variable_canon handles general complex variable", {
  z <- Variable(2, complex = TRUE)
  real2imag <- new.env(hash = TRUE, parent = emptyenv())
  imag_id <- 99999L
  assign(as.character(z@id), imag_id, envir = real2imag)

  result <- c2r_variable_canon(z, list(), list(), real2imag)
  expect_false(is.null(result[[1L]]))
  expect_false(is.null(result[[2L]]))
  ## Real part has same var_id as original
  expect_equal(result[[1L]]@id, z@id)
  ## Imag part has the mapped var_id
  expect_equal(result[[2L]]@id, imag_id)
})

## @cvxpy NONE
test_that("c2r_variable_canon handles imaginary variable", {
  z <- Variable(2, imag = TRUE)
  real2imag <- new.env(hash = TRUE, parent = emptyenv())
  imag_id <- 88888L
  assign(as.character(z@id), imag_id, envir = real2imag)

  result <- c2r_variable_canon(z, list(), list(), real2imag)
  expect_null(result[[1L]])
  expect_false(is.null(result[[2L]]))
  expect_equal(result[[2L]]@id, imag_id)
})

## @cvxpy NONE
test_that("c2r_variable_canon handles Hermitian variable", {
  n <- 3L
  Z <- Variable(c(n, n), hermitian = TRUE)
  real2imag <- new.env(hash = TRUE, parent = emptyenv())
  imag_id <- 77777L
  assign(as.character(Z@id), imag_id, envir = real2imag)

  result <- c2r_variable_canon(Z, list(), list(), real2imag)
  ## Real part: (n, n) symmetric variable
  expect_false(is.null(result[[1L]]))
  expect_true(S7_inherits(result[[1L]], Variable))
  expect_equal(result[[1L]]@shape, c(n, n))
  expect_true(isTRUE(result[[1L]]@attributes$symmetric))

  ## Imag part: skew-symmetric matrix built from n*(n-1)/2 variables
  expect_false(is.null(result[[2L]]))
})

# ── Affine canons ─────────────────────────────────────────────────

## @cvxpy NONE
test_that("c2r_separable_canon with all real args", {
  x <- Variable(2)
  expr <- -x  # NegExpression
  result <- c2r_separable_canon(expr, list(x), list(NULL), NULL)
  expect_false(is.null(result[[1L]]))
  expect_null(result[[2L]])
})

## @cvxpy NONE
test_that("c2r_separable_canon with mixed real/imag args", {
  x <- Variable(2)
  y <- Variable(2)
  expr <- x + y  # AddExpression
  result <- c2r_separable_canon(expr, list(x, NULL), list(NULL, y), NULL)
  expect_false(is.null(result[[1L]]))
  expect_false(is.null(result[[2L]]))
})

## @cvxpy NONE
test_that("c2r_conj_canon negates imaginary part", {
  x <- Variable(2)
  y <- Variable(2)
  result <- c2r_conj_canon(NULL, list(x), list(y), NULL)
  expect_identical(result[[1L]], x)
  ## Imaginary part is negated
  expect_true(S7_inherits(result[[2L]], NegExpression))
})

## @cvxpy NONE
test_that("c2r_real_canon returns real arg", {
  x <- Variable(2)
  result <- c2r_real_canon(NULL, list(x), list(NULL), NULL)
  expect_identical(result[[1L]], x)
  expect_null(result[[2L]])
})

## @cvxpy NONE
test_that("c2r_imag_canon returns imag arg as real", {
  y <- Variable(2)
  result <- c2r_imag_canon(NULL, list(NULL), list(y), NULL)
  expect_identical(result[[1L]], y)
  expect_null(result[[2L]])
})

# ── Binary canon ──────────────────────────────────────────────────

## @cvxpy NONE
test_that("c2r_binary_canon handles (a+bi)(c+di)", {
  a <- Constant(matrix(c(1, 0), 2, 1))
  b <- Constant(matrix(c(0, 1), 2, 1))
  c_ <- Constant(matrix(c(2, 0), 2, 1))
  d_ <- Constant(matrix(c(0, 3), 2, 1))

  ## Multiply atom — need some expression with 2 args
  expr <- a * c_  # just for shape
  result <- c2r_binary_canon(expr, list(a, c_), list(b, d_), NULL)
  ## real: a*c - b*d; imag: a*d + b*c
  expect_false(is.null(result[[1L]]))
  expect_false(is.null(result[[2L]]))
})

# ── Constraint canons ─────────────────────────────────────────────

## @cvxpy NONE
test_that("c2r_equality_canon with real-only args", {
  x <- Variable(2)
  y <- Variable(2)
  con <- Equality(x, y)
  result <- c2r_equality_canon(con, list(x, y), list(NULL, NULL), NULL)
  expect_length(result[[1L]], 1)
  expect_null(result[[2L]])
})

## @cvxpy NONE
test_that("c2r_equality_canon with complex args", {
  x_r <- Variable(2)
  y_r <- Variable(2)
  x_i <- Variable(2)
  y_i <- Variable(2)
  con <- Equality(Variable(2, complex = TRUE), Variable(2, complex = TRUE))
  real2imag <- new.env(hash = TRUE, parent = emptyenv())
  assign(as.character(con@id), 55555L, envir = real2imag)

  result <- c2r_equality_canon(con, list(x_r, y_r), list(x_i, y_i), real2imag)
  expect_length(result[[1L]], 1)  # real equality
  expect_length(result[[2L]], 1)  # imag equality
})

## @cvxpy NONE
test_that("c2r_psd_canon with imag part creates block matrix", {
  re <- Variable(c(2L, 2L))
  im <- Variable(c(2L, 2L))
  con <- PSD(Variable(c(2L, 2L), complex = TRUE))
  result <- c2r_psd_canon(con, list(re), list(im), NULL)
  ## Returns list of constraints
  expect_true(is.list(result[[1L]]))
  expect_length(result[[1L]], 1)
  ## The new PSD constraint should be 4x4
  new_con <- result[[1L]][[1L]]
  expect_equal(new_con@shape, c(4L, 4L))
})

# ── Full reduction apply ─────────────────────────────────────────

## @cvxpy NONE
test_that("Complex2Real apply produces real-only problem", {
  z <- Variable(2, complex = TRUE)
  prob <- Problem(Minimize(sum_entries(Re(z))), list(Im(z) == 1))

  c2r <- Complex2Real()
  result <- reduction_apply(c2r, prob)
  new_prob <- result[[1L]]
  inv_data <- result[[2L]]

  ## New problem should not contain any complex leaves
  expect_false(complex2real_accepts(new_prob))
  ## Inverse data should have real2imag mapping
  expect_true(!is.null(inv_data$real2imag))
})

## @cvxpy NONE
test_that("Complex2Real apply with Hermitian variable", {
  n <- 2L
  Z <- Variable(c(n, n), hermitian = TRUE)
  prob <- Problem(Minimize(Trace(Re(Z))), list(Z %>>% 0))

  c2r <- Complex2Real()
  result <- reduction_apply(c2r, prob)
  new_prob <- result[[1L]]

  ## New problem should not contain any complex leaves
  expect_false(complex2real_accepts(new_prob))
})

# ── Abs canon ─────────────────────────────────────────────────────

## @cvxpy NONE
test_that("c2r_abs_canon handles purely real", {
  x <- Variable(c(2L, 1L))
  result <- c2r_abs_canon(NULL, list(x), list(NULL), NULL)
  expect_false(is.null(result[[1L]]))
  expect_null(result[[2L]])
})

## @cvxpy NONE
test_that("c2r_abs_canon handles complex", {
  x_r <- Variable(c(2L, 1L))
  x_i <- Variable(c(2L, 1L))
  result <- c2r_abs_canon(NULL, list(x_r), list(x_i), NULL)
  expect_false(is.null(result[[1L]]))
  expect_null(result[[2L]])
  ## Result should have same shape as input
  expect_equal(result[[1L]]@shape, x_r@shape)
})

# ── Registry ──────────────────────────────────────────────────────

## @cvxpy NONE
test_that("C2R canonicalization S7 methods are registered", {
  expect_true(has_c2r_canon(AddExpression(list(Variable(), Variable()))))
  expect_true(has_c2r_canon(Variable()))
  expect_true(has_c2r_canon(Constant(1)))
  expect_true(has_c2r_canon(Equality(Variable(), Constant(1))))
  expect_true(has_c2r_canon(Abs(Variable())))
  expect_true(has_c2r_canon(QuadForm(Variable(2), Constant(diag(2)))))
  ## Default returns FALSE for unregistered types
  expect_false(has_c2r_canon(Minimize(Variable())))
})
