context("test-g01-grad")
a <- Variable(name = "a")

x <- Variable(2, name = "x")
y <- Variable(2, name = "y")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("Test gradient for affine_prod", {
  skip_on_cran()
  value(C) <- rbind(c(1,-2), c(3,4), c(-1,-3))
  value(A) <- rbind(c(3,2), c(-5,1))
  expect_warning(expr <- C %*% A)

  expect_equivalent(as.matrix(grad(expr)[[as.character(C@id)]]), rbind(c(3,0,0,2,0,0), c(0,3,0,0,2,0), c(0,0,3,0,0,2),
                                                    c(-5,0,0,1,0,0), c(0,-5,0,0,1,0), c(0,0,-5,0,0,1)))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), rbind(c(1,3,-1,0,0,0), c(-2,4,-3,0,0,0), c(0,0,0,1,3,-1),
                                                   c(0,0,0,-2,4,-3)))
})

test_that("Test gradient for p_norm", {
  skip_on_cran()
  value(x) <- c(-1,0)
  expr <- p_norm(x, 1)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(-1,0)))

  value(x) <- c(0,10)
  expr <- p_norm(x, 1)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(0,1)))

  value(x) <- c(-3,4)
  expr <- p_norm(x, 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(-3.0/5,4.0/5)))

  value(x) <- c(-1,2)
  expr <- p_norm(x, 0.5)
  expect_true(is.na(grad(expr)[[as.character(x@id)]]))

  value(x) <- c(0,0)
  expr <- p_norm(x, 0.5)
  expect_true(is.na(grad(expr)[[as.character(x@id)]]))

  value(x) <- c(0,0)
  expr <- p_norm(x, 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(0,0)))

  value(A) <- rbind(c(2,-2), c(2,2))
  expr <- p_norm(A, 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(0.5,0.5,-0.5,0.5)))

  value(A) <- rbind(c(3,-3), c(4,4))
  expr <- p_norm(A, 2, axis = 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), rbind(c(0.6,0), c(0.8,0), c(0,-0.6), c(0,0.8)))

  value(A) <- rbind(c(3,-4), c(4,3))
  expr <- p_norm(A, 2, axis = 1)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), rbind(c(0.6,0), c(0,0.8), c(-0.8,0), c(0,0.6)))

  value(A) <- rbind(c(3,-4), c(4,3))
  expr <- p_norm(A, 0.5)
  expect_true(is.na(grad(expr)[[as.character(A@id)]]))
})

test_that("Test gradient for log_sum_exp", {
  skip_on_cran()
  value(x) <- c(0,1)
  expr <- log_sum_exp(x)
  e <- exp(1)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(1.0/(1+e), e/(1+e))))

  value(A) <- rbind(c(0,1), c(-1,0))
  expr <- log_sum_exp(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(1.0/(2+e+1.0/e), 1.0/e/(2+e+1.0/e), e/(2+e+1.0/e), 1.0/(2+e+1.0/e))))

  value(A) <- rbind(c(0,1), c(-1,0))
  expr <- log_sum_exp(A, axis = 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), cbind(c(1.0/(1+1.0/e), 1.0/e/(1+1.0/e), 0, 0), c(0, 0, e/(1+e), 1.0/(1+e))))
})

test_that("Test gradient for geo_mean", {
  skip_on_cran()
  value(x) <- c(1,2)
  expr <- geo_mean(x)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(sqrt(2)/2, 1.0/2/sqrt(2))))

  value(x) <- c(0,2)
  expr <- geo_mean(x)
  expect_true(is.na(grad(expr)[[as.character(x@id)]]))

  value(x) <- c(1,2)
  expr <- geo_mean(x, c(1,0))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(1,0)))

  # No exception for single weight
  value(x) <- c(-1,2)
  expr <- geo_mean(x, c(1,0))
  expect_true(is.na(grad(expr)[[as.character(x@id)]]))
})

test_that("Test gradient for lambda_max", {
  skip_on_cran()
  value(A) <- cbind(c(2,0), c(0,1))
  expr <- lambda_max(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(1,0,0,0)))

  value(A) <- cbind(c(1,0), c(0,2))
  expr <- lambda_max(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(0,0,0,1)))

  value(A) <- cbind(c(1,0), c(0,1))
  expr <- lambda_max(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(0,0,0,1)))
})

test_that("Test gradient for matrix_frac", {
  skip_on_cran()
  value(A) <- diag(2)
  value(B) <- diag(2)
  expr <- matrix_frac(A, B)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(2,0,0,2)))
  expect_equivalent(as.matrix(grad(expr)[[as.character(B@id)]]), matrix(c(-1,0,0,-1)))

  value(B) <- matrix(0, nrow = 2, ncol = 2)
  expr <- matrix_frac(A, B)
  expect_true(is.na(grad(expr)[[as.character(A@id)]]))
  expect_true(is.na(grad(expr)[[as.character(B@id)]]))

  value(x) <- c(2,3)
  value(A) <- diag(2)
  expr <- matrix_frac(x, A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(4,6)))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(-4,-6,-6,-9)))
})

test_that("Test gradient for norm_nuc", {
  skip_on_cran()
  value(A) <- cbind(c(10,4), c(4,30))
  expr <- norm_nuc(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(1,0,0,1)))
})

test_that("Test gradient for log_det", {
  skip_on_cran()
  value(A) <- 2*diag(rep(1,2))
  expr <- log_det(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), 1.0/2*diag(2))

  mat <- rbind(c(1,2), c(3,5))
  value(A) <- t(mat) %*% mat
  expr <- log_det(A)
  val <- t(base::solve(value(A)))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)

  value(A) <- matrix(0, nrow = 2, ncol = 2)
  expr <- log_det(A)
  expect_true(is.na(grad(expr)[[as.character(A@id)]]))

  value(A) <- -rbind(c(1,2), c(3,4))
  expr <- log_det(A)
  expect_true(is.na(grad(expr)[[as.character(A@id)]]))
})

test_that("Test gradient for quad_over_lin", {
  skip_on_cran()
  value(x) <- c(1,2)
  value(a) <- 2
  expr <- quad_over_lin(x, a)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(1,2)))
  expect_equal(grad(expr)[[as.character(a@id)]], -1.25)

  value(a) <- 0
  expr <- quad_over_lin(x, a)
  expect_true(is.na(grad(expr)[[as.character(x@id)]]))
  expect_true(is.na(grad(expr)[[as.character(a@id)]]))

  value(A) <- diag(rep(1,2))
  value(a) <- 2
  expr <- quad_over_lin(A, a)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(1,0,0,1)))
  expect_equal(grad(expr)[[as.character(a@id)]], -0.5)

  value(x) <- c(1,2)
  value(a) <- 2
  value(y) <- c(1,2)
  value(a) <- 2
  expr <- quad_over_lin(x, a) + quad_over_lin(y, a)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(1,2)))
  expect_equivalent(as.matrix(grad(expr)[[as.character(y@id)]]), matrix(c(1,2)))
  expect_equal(grad(expr)[[as.character(a@id)]], -2.5)
})

test_that("Test gradient for max_entries", {
  skip_on_cran()
  value(x) <- c(2,1)
  expr <- max_entries(x)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(c(1,0)))

  value(A) <- rbind(c(1,2), c(4,3))
  expr <- max_entries(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(0,1,0,0)))

  value(A) <- rbind(c(1,2), c(4,3))
  expr <- max_entries(A, axis = 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), rbind(c(0,0), c(1,0), c(0,0), c(0,1)))

  value(A) <- rbind(c(1,2), c(4,3))
  expr <- max_entries(A, axis = 1)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), rbind(c(0,0), c(0,1), c(1,0), c(0,0)))
})

test_that("Test sigma_max", {
  skip_on_cran()
  value(A) <- cbind(c(1,0), c(0,2))
  expr <- sigma_max(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(0,0,0,1)))

  value(A) <- cbind(c(1,0), c(0,1))
  expr <- sigma_max(A)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(1,0,0,0)))
})

test_that("Test sum_largest", {
  skip_on_cran()
  value(A) <- cbind(c(4,3), c(2,1))
  expr <- sum_largest(A, 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(1,0,1,0)))

  value(A) <- cbind(c(1,2), c(3,0.5))
  expr <- sum_largest(A, 2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), matrix(c(0,1,1,0)))
})

test_that("Test abs", {
  skip_on_cran()
  value(A) <- cbind(c(1,2), c(-1,0))
  expr <- abs(A)
  val <- diag(c(1,1,-1,0))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)
})

test_that("Test linearize method", {
  skip_on_cran()
  # Affine
  value(x) <- c(1,2)
  expr <- (2*x - 5)[1]
  lin_expr <- linearize(expr)
  value(x) <- c(55,22)
  expr <- (2*x - 5)[1]
  lin_expr <- linearize(expr)
  expect_equal(value(lin_expr), value(expr))
  value(x) <- c(-1,-5)
  expr <- (2*x - 5)[1]
  lin_expr <- linearize(expr)
  expect_equal(value(lin_expr), value(expr))

  # Convex
  expr <- A^2 + 5
  expect_error(linearize(expr))

  value(A) <- cbind(c(1,2), c(3,4))
  expr <- A^2 + 5
  lin_expr <- linearize(expr)
  manual <- value(expr) + 2*reshape_expr(value(diag(vec(A))) %*% vec(A - value(A)), c(2, 2))
  expect_equivalent(as.matrix(value(lin_expr)), value(expr))

  value(A) <- rbind(c(-5,-5), c(8.2,4.4))
  expr <- A^2 + 5
  lin_expr <- linearize(expr)
  manual <- value(expr) + 2*reshape_expr(value(diag(vec(A))) %*% vec(A - value(A)), c(2, 2))
  expect_true(all(value(lin_expr) <= value(expr)))
  expect_equivalent(as.matrix(value(lin_expr)), value(manual))

  # Concave
  value(x) <- c(1,2)
  expr <- log(x)/2
  lin_expr <- linearize(expr)
  manual <- value(expr) + value(diag(0.5*x^-1)) %*% (x - value(x))
  expect_equivalent(as.matrix(value(lin_expr)), value(expr))

  value(x) <- c(3,4.4)
  expr <- log(x)/2
  lin_expr <- linearize(expr)
  manual <- value(expr) + value(diag(0.5*x^-1)) %*% (x - value(x))
  expect_true(all(value(lin_expr) >= value(expr)))
  expect_equivalent(as.matrix(value(lin_expr)), value(manual))
})

test_that("Test gradient for log", {
  skip_on_cran()
  value(a) <- 2
  expr <- log(a)
  expect_equal(grad(expr)[[as.character(a@id)]], 1.0/2)

  value(a) <- 3
  expr <- log(a)
  expect_equal(grad(expr)[[as.character(a@id)]], 1.0/3)

  value(a) <- -1
  expr <- log(a)
  expect_true(is.na(grad(expr)[[as.character(a@id)]]))

  value(x) <- c(3,4)
  expr <- log(x)
  val <- diag(c(1/3,1/4))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)

  value(x) <- c(-1e-9,4)
  expr <- log(x)
  expect_true(is.na(grad(expr)[[as.character(x@id)]]))

  value(A) <- cbind(c(1,2), c(3,4))
  expr <- log(A)
  val <- diag(c(1, 1/2, 1/3, 1/4))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)
})

test_that("Test domain for log1p", {
  skip_on_cran()
  value(a) <- 2
  expr <- log1p(a)
  expect_equal(grad(expr)[[as.character(a@id)]], 1.0/3)

  value(a) <- 3
  expr <- log1p(a)
  expect_equal(grad(expr)[[as.character(a@id)]], 1.0/4)

  value(a) <- -1
  expr <- log1p(a)
  expect_true(is.na(grad(expr)[[as.character(a@id)]]))

  value(x) <- c(3,4)
  expr <- log1p(x)
  val <- diag(c(1/4,1/5))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)

  value(x) <- c(-1e-9-1,4)
  expr <- log1p(x)
  expect_true(is.na(grad(expr)[[as.character(x@id)]]))

  value(A) <- cbind(c(1,2), c(3,4))
  expr <- log1p(A)
  val <- diag(c(1/2, 1/3, 1/4, 1/5))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)
})

test_that("Test domain for entr", {
  skip_on_cran()
  value(a) <- 2
  expr <- entr(a)
  expect_equal(grad(expr)[[as.character(a@id)]], -log(2)-1)

  value(a) <- 3
  expr <- entr(a)
  expect_equal(grad(expr)[[as.character(a@id)]], -(log(3)+1))

  value(a) <- -1
  expr <- entr(a)
  expect_true(is.na(grad(expr)[[as.character(a@id)]]))

  value(x) <- c(3,4)
  expr <- entr(x)
  val <- diag(-(log(c(3,4)) + 1))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)

  value(x) <- c(-1e-9,4)
  expr <- entr(x)
  expect_true(is.na(grad(expr)[[as.character(x@id)]]))

  value(A) <- cbind(c(1,2), c(3,4))
  expr <- entr(A)
  val <- diag(-(log(1:4)+1))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)
})

test_that("Test domain for exp", {
  skip_on_cran()
  value(a) <- 2
  expr <- exp(a)
  expect_equal(grad(expr)[[as.character(a@id)]], exp(2))

  value(a) <- 3
  expr <- exp(a)
  expect_equal(grad(expr)[[as.character(a@id)]], exp(3))

  value(a) <- -1
  expr <- exp(a)
  expect_equal(grad(expr)[[as.character(a@id)]], exp(-1))

  value(x) <- c(3,4)
  expr <- exp(x)
  val <- diag(exp(c(3,4)))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)

  value(x) <- c(-1e-9,4)
  expr <- exp(x)
  val <- diag(exp(c(-1e-9,4)))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)

  value(A) <- cbind(c(1,2), c(3,4))
  expr <- exp(A)
  val <- diag(exp(1:4))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)
})

test_that("Test domain for logistic", {
  skip_on_cran()
  value(a) <- 2
  expr <- logistic(a)
  expect_equal(grad(expr)[[as.character(a@id)]], exp(2)/(1+exp(2)))

  value(a) <- 3
  expr <- logistic(a)
  expect_equal(grad(expr)[[as.character(a@id)]], exp(3)/(1+exp(3)))

  value(a) <- -1
  expr <- logistic(a)
  expect_equal(grad(expr)[[as.character(a@id)]], exp(-1)/(1+exp(-1)))

  value(x) <- c(3,4)
  expr <- logistic(x)
  val <- diag(exp(c(3,4))/(1+exp(c(3,4))))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)

  value(x) <- c(-1e-9,4)
  expr <- logistic(x)
  val <- diag(exp(c(-1e-9,4))/(1+exp(c(-1e-9,4))))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)

  value(A) <- cbind(c(1,2), c(3,4))
  expr <- logistic(A)
  val <- diag(exp(1:4)/(1+exp(1:4)))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)
})

test_that("Test domain for huber", {
  skip_on_cran()
  value(a) <- 2
  expr <- CVXR::huber(a)
  expect_equal(grad(expr)[[as.character(a@id)]], 2)

  value(a) <- 3
  expr <- CVXR::huber(a, M = 2)
  expect_equal(grad(expr)[[as.character(a@id)]], 4)

  value(a) <- -1
  expr <- CVXR::huber(a, M = 2)
  expect_equal(grad(expr)[[as.character(a@id)]], -2)

  value(x) <- c(3,4)
  expr <- CVXR::huber(x)
  val <- diag(c(2,2))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)

  value(x) <- c(-1e-9,4)
  expr <- CVXR::huber(x)
  val <- diag(c(0,2))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)

  value(A) <- cbind(c(1,2), c(3,4))
  expr <- CVXR::huber(A, M = 3)
  val <- diag(c(2,4,6,6))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)
})

test_that("Test domain for kl_div", {
  skip_on_cran()
  b <- Variable()
  value(a) <- 2
  value(b) <- 4
  expr <- kl_div(a, b)
  expect_equal(grad(expr)[[as.character(a@id)]], log(2/4))
  expect_equal(grad(expr)[[as.character(b@id)]], 1-(2/4))

  value(a) <- 3
  value(b) <- 0
  expr <- kl_div(a, b)
  expect_true(is.na(grad(expr)[[as.character(a@id)]]))
  expect_true(is.na(grad(expr)[[as.character(b@id)]]))

  value(a) <- -1
  value(b) <- 2
  expr <- kl_div(a, b)
  expect_true(is.na(grad(expr)[[as.character(a@id)]]))
  expect_true(is.na(grad(expr)[[as.character(b@id)]]))

  y <- Variable(2)
  value(x) <- c(3,4)
  value(y) <- c(5,8)
  expr <- kl_div(x, y)
  val <- diag(log(c(3,4)) - log(c(5,8)))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  val <- diag(c(1-3/5,1-4/8))
  expect_equivalent(as.matrix(grad(expr)[[as.character(y@id)]]), val)

  value(x) <- c(-1e-9,4)
  value(y) <- c(1,2)
  expr <- kl_div(x, y)
  expect_true(is.na(grad(expr)[[as.character(x@id)]]))
  expect_true(is.na(grad(expr)[[as.character(y@id)]]))

  value(A) <- cbind(c(1,2), c(3,4))
  value(B) <- cbind(c(5,1), c(3.5,2.3))
  expr <- kl_div(A, B)
  div <- as.numeric(value(A) / value(B))
  val <- diag(log(div))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)
  val <- diag(1-div)
  expect_equivalent(as.matrix(grad(expr)[[as.character(B@id)]]), val)
})

test_that("Test domain for max_elemwise", {
  skip_on_cran()
  b <- Variable()
  value(a) <- 2
  value(b) <- 4
  expr <- max_elemwise(a, b)
  expect_equal(grad(expr)[[as.character(a@id)]], 0)
  expect_equal(grad(expr)[[as.character(b@id)]], 1)

  value(a) <- 3
  value(b) <- 0
  expr <- max_elemwise(a, b)
  expect_equal(grad(expr)[[as.character(a@id)]], 1)
  expect_equal(grad(expr)[[as.character(b@id)]], 0)

  value(a) <- -1
  value(b) <- 2
  expr <- max_elemwise(a, b)
  expect_equal(grad(expr)[[as.character(a@id)]], 0)
  expect_equal(grad(expr)[[as.character(b@id)]], 1)

  y <- Variable(2)
  value(x) <- c(3,4)
  value(y) <- c(5,-5)
  expr <- max_elemwise(x, y)
  val <- diag(c(0,1))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  val <- diag(c(1,0))
  expect_equivalent(as.matrix(grad(expr)[[as.character(y@id)]]), val)

  value(x) <- c(-1e-9,4)
  value(y) <- c(1,4)
  expr <- max_elemwise(x, y)
  val <- diag(c(0,1))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  val <- diag(c(1,0))
  expect_equivalent(as.matrix(grad(expr)[[as.character(y@id)]]), val)

  value(A) <- cbind(c(1,2), c(3,4))
  value(B) <- cbind(c(5,1), c(3,2.3))
  expr <- max_elemwise(A, B)
  div <- as.numeric(value(A) / value(B))
  val <- diag(c(0,1,1,1))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)
  val <- diag(c(1,0,0,0))
  expect_equivalent(as.matrix(grad(expr)[[as.character(B@id)]]), val)
})

test_that("Test domain for min_elemwise", {
  skip_on_cran()
  b <- Variable()
  value(a) <- 2
  value(b) <- 4
  expr <- min_elemwise(a, b)
  expect_equal(grad(expr)[[as.character(a@id)]], 1)
  expect_equal(grad(expr)[[as.character(b@id)]], 0)

  value(a) <- 3
  value(b) <- 0
  expr <- min_elemwise(a, b)
  expect_equal(grad(expr)[[as.character(a@id)]], 0)
  expect_equal(grad(expr)[[as.character(b@id)]], 1)

  value(a) <- -1
  value(b) <- 2
  expr <- min_elemwise(a, b)
  expect_equal(grad(expr)[[as.character(a@id)]], 1)
  expect_equal(grad(expr)[[as.character(b@id)]], 0)

  y <- Variable(2)
  value(x) <- c(3,4)
  value(y) <- c(5,-5)
  expr <- min_elemwise(x, y)
  val <- diag(c(1,0))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  val <- diag(c(0,1))
  expect_equivalent(as.matrix(grad(expr)[[as.character(y@id)]]), val)

  value(x) <- c(-1e-9,4)
  value(y) <- c(1,4)
  expr <- min_elemwise(x, y)
  val <- diag(c(1,1))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)
  val <- diag(c(0,0))
  expect_equivalent(as.matrix(grad(expr)[[as.character(y@id)]]), val)

  value(A) <- cbind(c(1,2), c(3,4))
  value(B) <- cbind(c(5,1), c(3,2.3))
  expr <- min_elemwise(A, B)
  div <- as.numeric(value(A) / value(B))
  val <- diag(c(1,0,1,0))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)
  val <- diag(c(0,1,0,1))
  expect_equivalent(as.matrix(grad(expr)[[as.character(B@id)]]), val)
})

test_that("Test domain for power", {
  skip_on_cran()
  value(a) <- 2
  expr <- sqrt(a)
  expect_equal(grad(expr)[[as.character(a@id)]], 0.5/sqrt(2))

  value(a) <- 3
  expr <- sqrt(a)
  expect_equal(grad(expr)[[as.character(a@id)]], 0.5/sqrt(3))

  value(a) <- -1
  expr <- sqrt(a)
  expect_true(is.na(grad(expr)[[as.character(a@id)]]))

  value(x) <- c(3,4)
  expr <- x^3
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), rbind(c(27,0), c(0,48)))

  value(x) <- c(-1e-9,4)
  expr <- x^3
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), rbind(c(0,0), c(0,48)))

  value(A) <- cbind(c(1,-2), c(3,4))
  expr <- A^2
  val <- diag(c(2,-4,6,8))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)

  # Constant
  expr <- a^0
  expect_equal(grad(expr)[[as.character(a@id)]], 0)

  expr <- x^0
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), matrix(0, nrow = 2, ncol = 2))
})

# test_that("Test grad for partial minimization/maximization problems", {
  # for(obj in list(Minimize(a^-1), Maximize(Entr(a)))) {
  #   prob <- Problem(obj, list(x + a >= c(5,8)))
  #
  #   # Optimize over nothing
  #   expr <- partial_optimize(prob, dont_opt_vars = list(x, a))
  #   value(a) <- NA
  #   value(x) <- NA
  #   grad <- grad(expr)
  #   expect_true(is.na(grad[[as.character(a@id)]]))
  #   expect_true(is.na(grad[[as.character(x@id)]]))
  #
  #   # Outside domain
  #   value(a) <- 1.0
  #   value(x) <- c(5,5)
  #   grad <- grad(expr)
  #   expect_true(is.na(grad[[as.character(a@id)]]))
  #   expect_true(is.na(grad[[as.character(x@id)]]))
  #
  #   value(a) <- 1
  #   value(x) <- c(10,10)
  #   grad <- grad(expr)
  #   expect_equal(grad[[as.character(a@id)]], grad(obj@.args[[1]])[[as.character(a@id)]])
  #   expect_equal(as.matrix(grad[[as.character(x@id)]]), rep(0,4))
  #
  #   # Optimize over x
  #   expr <- partial_optimize(prob, opt_vars = list(x))
  #   value(a) <- 1
  #   grad <- grad(expr)
  #   expect_equal(grad[[as.character(a@id)]], grad(obj@.args[[1]])[[as.character(a@id)]] + 0)
  #
  #   # Optimize over a
  #   fix_prob <- Problem(obj, list(x + a >= c(5,8), x == 0))
  #   result <- solve(fix_prob)
  #   dual_val <- fix_prob@constraints[1]@dual_variable@value
  #   expr <- partial_optimize(prob, opt_vars = list(a))
  #   value(x) <- c(0,0)
  #   grad <- grad(expr)
  #   expect_equal(as.matrix(grad[[as.character(x@id)]]), dual_val)
  #
  #   # Optimize over x and a
  #   expr <- partial_optimize(prob, opt_vars = list(x, a))
  #   grad <- grad(expr)
  #   expect_equal(grad, list())
  # }
# })

test_that("Test grad for affine atoms", {
  skip_on_cran()
  value(a) <- 2
  expr <- -a
  expect_equal(grad(expr)[[as.character(a@id)]], -1)

  value(a) <- 2
  expr <- 2*a
  expect_equal(grad(expr)[[as.character(a@id)]], 2)

  value(a) <- 2
  expr <- a/2
  expect_equal(grad(expr)[[as.character(a@id)]], 0.5)

  value(x) <- c(3,4)
  expr <- -x
  val <- -diag(c(1,1))
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)

  value(A) <- cbind(c(1,2), c(3,4))
  expr <- -A
  val <- -diag(rep(1,4))
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)

  value(A) <- cbind(c(1,2), c(3,4))
  expr <- A[1,2]
  val <- matrix(0, nrow = 4, ncol = 1)
  val[3] <- 1
  expect_equivalent(as.matrix(grad(expr)[[as.character(A@id)]]), val)

  z <- Variable(3)
  value(x) <- c(1,2)
  value(z) <- c(1,2,3)
  expr <- vstack(x, z)
  val <- matrix(0, nrow = 2, ncol = 5)
  val[,1:2] <- diag(2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)

  val <- matrix(0, nrow = 3, ncol = 5)
  val[,3:ncol(val)] <- diag(3)
  expect_equivalent(as.matrix(grad(expr)[[as.character(z@id)]]), val)

  # Cumulative sum
  value(x) <- c(1,2)
  expr <- cumsum_axis(x)
  val <- matrix(1, nrow = 2, ncol = 2)
  val[2,1] <- 0
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)

  value(x) <- c(1,2)
  expr <- cumsum_axis(x, axis = 1)
  val <- diag(2)
  expect_equivalent(as.matrix(grad(expr)[[as.character(x@id)]]), val)
})
