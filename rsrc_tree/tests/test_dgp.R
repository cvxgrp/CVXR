context("test_dgp")
TOL <- 1e-6

test_that("test product", {
  skip_on_cran()
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prod <- x*y
  expect_true(is_dgp(prod))
  expect_true(is_log_log_convex(prod))
  expect_true(is_log_log_concave(prod))

  prod <- prod*prod
  expect_true(is_dgp(prod))
  expect_true(is_log_log_convex(prod))
  expect_true(is_log_log_concave(prod))

  prod <- 5.0*prod
  expect_true(is_dgp(prod))
  expect_true(is_log_log_convex(prod))
  expect_true(is_log_log_concave(prod))

  prod <- -5.0*prod
  expect_false(is_dgp(prod))
  expect_false(is_log_log_convex(prod))
  expect_false(is_log_log_concave(prod))
})

test_that("test product with unconstrained variables is not DGP", {
  skip_on_cran()
  x <- Variable()
  y <- Variable()
  prod <- x*y
  expect_false(is_dgp(prod))
  expect_false(is_log_log_convex(prod))
  expect_false(is_log_log_concave(prod))

  z <- Variable(pos = TRUE)
  prod <- x*z
  expect_false(is_dgp(prod))
  expect_false(is_log_log_convex(prod))
  expect_false(is_log_log_concave(prod))
})

test_that("test division", {
  skip_on_cran()
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  div <- x/y

  expect_true(is_log_log_affine(div))

  posynomial <- 5.0*x*y + 1.2*y*y
  div <- posynomial/(3.0*x*y^(-0.1))
  expect_true(is_log_log_convex(div))
  expect_false(is_log_log_concave(div))
  expect_true(is_dgp(div))

  div <- posynomial/(3.0*x + y)
  expect_false(is_log_log_convex(div))
  expect_false(is_log_log_concave(div))
  expect_false(is_dgp(div))
})

test_that("test add", {
  skip_on_cran()
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  expr <- x + y
  expect_true(is_dgp(expr))
  expect_true(is_log_log_convex(expr))
  expect_false(is_log_log_concave(expr))

  posynomial <- 5.0*x*y + 1.2*y*y
  expect_true(is_dgp(posynomial))
  expect_true(is_log_log_convex(posynomial))
})

test_that("test add with unconstrained variables is not DGP", {
  skip_on_cran()
  x <- Variable()
  y <- Variable(pos = TRUE)
  expr <- x + y
  expect_false(is_dgp(expr))
  expect_false(is_log_log_convex(expr))
  expect_false(is_log_log_concave(expr))

  posynomial <- 5.0*x*y + 1.2*y*y
  expect_false(is_dgp(posynomial))
  expect_false(is_log_log_convex(posynomial))
  expect_false(is_log_log_concave(posynomial))
})

test_that("test monomials", {
  skip_on_cran()
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  monomial <- 5.0*(x^(0.1)) * y^(-0.1) * z^3
  expect_true(is_dgp(monomial))
  expect_true(is_log_log_convex(monomial))
  expect_true(is_log_log_concave(monomial))

  monomial <- -1.0*monomial
  expect_false(is_dgp(monomial))
  expect_false(is_log_log_convex(monomial))
  expect_false(is_log_log_concave(monomial))
})

test_that("test elementwise maximum", {
  skip_on_cran()
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  monomial <- 5.0*(x^(0.1)) * y^(-0.1) * z^3
  posynomial <- 5.0*x*y + 1.2*y*y
  another_posynomial <- posynomial*posynomial
  expr <- max_elemwise(monomial, posynomial, another_posynomial)
  expect_true(is_dgp(expr))
  expect_true(is_log_log_convex(expr))
  expect_false(is_log_log_concave(expr))

  expr <- posynomial*expr
  expect_true(is_dgp(expr))
  expect_true(is_log_log_convex(expr))
  expect_false(is_log_log_concave(expr))

  expr <- posynomial*expr + expr
  expect_true(is_dgp(expr))
  expect_true(is_log_log_convex(expr))
})

test_that("test elementwise minimum", {
  skip_on_cran()
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  monomial <- 5.0*(x^(0.1)) * y^(-0.1) * z^3
  posynomial <- 5.0*x*y + 1.2*y*y
  another_posynomial <- posynomial*posynomial
  expr <- min_elemwise(monomial, 1/posynomial, 1/another_posynomial)
  expect_true(is_dgp(expr))
  expect_false(is_log_log_convex(expr))
  expect_true(is_log_log_concave(expr))

  expr <- (1/posynomial)*expr
  expect_true(is_dgp(expr))
  expect_false(is_log_log_convex(expr))
  expect_true(is_log_log_concave(expr))

  expr <- expr^2
  expect_true(is_dgp(expr))
  expect_false(is_log_log_convex(expr))
  expect_true(is_log_log_concave(expr))
})

test_that("test constant", {
  skip_on_cran()
  x <- Constant(1.0)
  expect_true(is_dgp(x))
  expect_false(is_dgp(-1.0*x))
})

test_that("test geometric mean", {
  skip_on_cran()
  x <- Variable(3, pos = TRUE)
  p <- c(1, 2, 0.5)
  expr <- geo_mean(x, p)
  expect_true(is_dgp(expr))
  expect_true(is_log_log_affine(expr))
  expect_true(is_log_log_convex(expr))
  expect_true(is_log_log_concave(expr))
})

test_that("test builtin sum", {
  x <- Variable(2, pos = TRUE)
  expect_true(is_log_log_convex(sum(x)))
})

test_that("test gmatmul", {
  x <- Variable(2, pos = TRUE)
  A <- Variable(2, 2)
  expect_error(gmatmul(A, x), "gmatmul(A, X) requires that A be constant", fixed = TRUE)
  
  x <- Variable(2)
  A <- matrix(1, nrow = 4, ncol = 2)
  expect_error(gmatmul(A, x), "gmatmul(A, X) requires that X be positive", fixed = TRUE)
  
  x <- Variable(3, pos = TRUE)
  A <- matrix(1, nrow = 4, ncol = 3)
  expr <- gmatmul(A, x)
  expect_true(is_dgp(expr))
  expect_true(is_log_log_affine(expr))
  expect_true(is_log_log_convex(expr))
  expect_true(is_log_log_concave(expr))
  expect_true(is_nonneg(expr))
  expect_true(is_incr(expr, 1))
  expect_true(is_decr(gmatmul(-A, x), 1))
  
  x <- Variable(2, 3, pos = TRUE)
  A <- rbind(c(2, -1), c(0, 3))
  expr <- gmatmul(A, x)
  expect_true(is_dgp(expr))
  expect_true(is_log_log_affine(expr))
  expect_true(is_log_log_convex(expr))
  expect_true(is_log_log_concave(expr))
  expect_false(is_incr(expr, 1))
  expect_false(is_decr(expr, 1))
})

test_that("test power sign", {
  x <- Variable(pos = TRUE)
  expect_true(is_nonneg(x^1))
  expect_false(is_nonpos(x^1))
})

test_that("test sparse constant not allowed", {
  sparse_matrix <- Constant(Matrix(c(1.0, 2.0), sparse = TRUE))
  expect_false(is_log_log_constant(sparse_matrix))
})
