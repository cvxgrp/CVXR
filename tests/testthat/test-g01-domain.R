context("test-g01-domain")
TOL <- 1e-6

a <- Variable(name = "a")

x <- Variable(2, name = "x")
y <- Variable(2, name = "y")
z <- Variable(3, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

# test_that("Test domain for partial minimization/maximization problems", {
#   for(obj in list(Minimize(a^-1), Maximize(log(a)))) {
#     prob <- Problem(obj, list(x + a >= c(5,8)))

    # Optimize over nothing
    # expr <- partial_optimize(prob, dont_opt_vars = list(x, a))
    # dom <- domain(expr)
    # constr <- list(a >= -100, x >= 0)
    # prob <- Problem(Minimize(sum(x + a)), c(dom, constr))
    # result <- solve(prob)
    # expect_equal(result$value, 13, tolerance = TOL)
    # expect_true(result$getValue(a) >= 0)
    # expect_true(all(value(x + a - c(5, 8)) >= -1e-3))

    # Optimize over x
    # expr <- partial_optimize(prob, opt_vars = list(x))
    # dom <- domain(expr)
    # constr <- list(a >= -100, x >= 0)
    # prob <- Problem(Minimize(sum(x + a)), c(dom, constr))
    # result <- solve(prob)
    # expect_equal(result$value, 0, tolerance = TOL)
    # expect_true(result$getValue(a) >= 0)
    # expect_equal(result$getValue(x), matrix(c(0,0)))

    # Optimize over x and a
    # expr <- partial_optimize(prob, opt_vars = list(x, a))
    # dom <- domain(expr)
    # constr <- list(a >= -100, x >= 0)
    # prob <- Problem(Minimize(sum(x + a)), c(dom, constr))
    # result <- solve(prob)
    # expect_equal(result$getValue(a), -100, tolerance = TOL)
    # expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = TOL)
#   }
# })

test_that("Test domain for geo_mean", {
  skip_on_cran()
  dom <- domain(geo_mean(x))
  prob <- Problem(Minimize(sum(x)), dom)
  result <- solve(prob)
  expect_equal(result$value, 0, tolerance = TOL)

  # No special case for only one weight
  dom <- domain(geo_mean(x, c(0,2)))
  dom <- c(dom, x >= -1)
  prob <- Problem(Minimize(sum(x)), dom)
  result <- solve(prob)
  expect_equal(result$getValue(x), matrix(c(-1,0)), tolerance = TOL)

  dom <- domain(geo_mean(z, c(0,1,1)))
  dom <- c(dom, z >= -1)
  prob <- Problem(Minimize(sum(z)), dom)
  result <- solve(prob)
  expect_equal(result$getValue(z), matrix(c(-1,0,0)))
})

test_that("Test domain for quad_over_lin", {
  skip_on_cran()
  dom <- domain(quad_over_lin(x, a))
  result <- solve(Problem(Minimize(a), dom))
  expect_equal(result$getValue(a), 0, tolerance = TOL)
})

test_that("Test domain for lambda_max", {
  skip_on_cran()
  dom <- domain(lambda_max(A))
  A0 <- rbind(c(1,2), c(3,4))
  result <- solve(Problem(Minimize(norm2(A-A0)), dom))
  expect_equal(result$getValue(A), rbind(c(1,2.5), c(2.5,4)), tolerance = TOL)
})

test_that("Test domain for p_norm", {
  skip_on_cran()
  dom <- domain(p_norm(a, -0.5))
  prob <- Problem(Minimize(a), dom)
  result <- solve(prob)
  expect_equal(result$value, 0, tolerance = TOL)
})

test_that("Test domain for log", {
  skip_on_cran()
  dom  <- domain(log(a))
  result <- solve(Problem(Minimize(a), dom))
  expect_equal(result$getValue(a), 0, tolerance = TOL)
})

test_that("Test domain for log1p", {
  skip_on_cran()
  dom <- domain(log1p(a))
  result <- solve(Problem(Minimize(a), dom))
  expect_equal(result$getValue(a), -1, tolerance = TOL)
})

test_that("Test domain for entr", {
  skip_on_cran()
  dom <- domain(entr(a))
  result <- solve(Problem(Minimize(a), dom))
  expect_equal(result$getValue(a), 0, tolerance = TOL)
})

test_that("Test domain for kl_div", {
  skip_on_cran()
  b <- Variable()
  dom <- domain(kl_div(a, b))
  result <- solve(Problem(Minimize(a + b), dom))
  expect_equal(result$getValue(a), 0, tolerance = TOL)
  expect_equal(result$getValue(b), 0, tolerance = TOL)
})

test_that("Test domain for power", {
  skip_on_cran()
  dom <- domain(sqrt(a))
  result <- solve(Problem(Minimize(a), dom))
  expect_equal(result$getValue(a), 0, tolerance = TOL)

  dom <- domain(a^2)
  result <- solve(Problem(Minimize(a), c(dom, a >= -100)))
  expect_equal(result$getValue(a), -100, tolerance = TOL)

  dom <- domain(a^-1)
  result <- solve(Problem(Minimize(a), c(dom, a >= -100)))
  expect_equal(result$getValue(a), 0, tolerance = TOL)

  dom <- domain(a^3)
  result <- solve(Problem(Minimize(a), c(dom, a >= -100)))
  expect_equal(result$getValue(a), 0, tolerance = TOL)
})

test_that("Test domain for log_det", {
  skip_on_cran()
  dom <- domain(log_det(A + diag(rep(1,2))))
  prob <- Problem(Minimize(sum(diag(A))), dom)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, -2, tolerance = 1e-3)
})

test_that("Test domain for matrix_frac", {
  skip_on_cran()
  dom <- domain(matrix_frac(x, A + diag(rep(1,2))))
  prob <- Problem(Minimize(sum(diag(A))), dom)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, -2, tolerance = 1e-3)
})
