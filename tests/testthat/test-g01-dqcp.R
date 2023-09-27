context("test-g01-dqcp")
TOL <- 1e-6
SOLVER <- "ECOS"

test_that("test basic with interval", {
  x <- Variable()
  expr <- ceil(x)
  
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconvex(expr))
  expect_true(is_quasiconcave(expr))
  expect_false(is_convex(expr))
  expect_false(is_concave(expr))
  expect_false(is_dcp(expr))
  expect_false(is_dgp(expr))
  
  problem <- Problem(Minimize(expr), list(x >= 12, x <= 17))
  expect_true(is_dqcp(problem))
  expect_false(is_dcp(problem))
  expect_false(is_dgp(problem))
  
  red <- Dqcp2Dcp(problem)
  reduced <- reduce(red)
  expect_true(is_dcp(reduced))
  expect_equal(length(parameters(reduced)), 1)
  soln <- bisect(reduced, low = 12, high = 17, solver = "SCS")
  expect_equal(soln@opt_val, 12.0, tolerance = 1e-3)
  
  result <- unpack_problem(problem, soln)
  expect_equal(soln@opt_val, result$value)
  expect_equal(result$getValue(x), 12.0, tolerance = 1e-3)
})

test_that("test basic without interval", {
  x <- Variable()
  expr <- ceil(x)
  
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconvex(expr))
  expect_true(is_quasiconcave(expr))
  expect_false(is_convex(expr))
  expect_false(is_concave(expr))
  expect_false(is_dcp(expr))
  expect_false(is_dgp(expr))
  
  problem <- Problem(Minimize(expr), list(x >= 12, x <= 17))
  expect_true(is_dqcp(problem))
  expect_false(is_dcp(problem))
  expect_false(is_dgp(problem))
  
  red <- Dqcp2Dcp(problem)
  reduced <- reduce(red)
  expect_true(is_dcp(reduced))
  expect_equal(length(parameters(reduced)), 1)
  soln <- bisect(reduced, solver = "SCS")
  expect_equal(soln@opt_val, 12.0, tolerance = 1e-3)
  
  result <- unpack_problem(problem, soln)
  expect_equal(soln@opt_val, result$value)
  expect_equal(result$getValue(x), 12.0, tolerance = 1e-3)
})

test_that("test basic solve", {
  x <- Variable()
  expr <- ceil(x)
  
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconvex(expr))
  expect_true(is_quasiconcave(expr))
  expect_false(is_convex(expr))
  expect_false(is_concave(expr))
  expect_false(is_dcp(expr))
  expect_false(is_dgp(expr))
  
  problem <- Problem(Minimize(expr), list(x >= 12, x <= 17))
  expect_true(is_dqcp(problem))
  expect_false(is_dcp(problem))
  expect_false(is_dgp(problem))
  result <- solve(problem, SOLVER, qcp = TRUE, low = 12, high = 17)
  expect_equal(result$value, 12.0, tolerance = 1e-3)
  expect_equal(result$getValue(x), 12.0, tolerance = 1e-3)
  
  problem <- clear_solution(problem)
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$value, 12.0, tolernace = 1e-3)
  expect_equal(result$getValue(x), 12.0, tolerance = 1e-3)
  
  problem <- clear_solution(problem)
  result <- solve(problem, SOLVER, qcp = TRUE, high = 17)
  expect_equal(result$value, 12.0, tolerance = 1e-3)
  expect_equal(result$getValue(x), 12.0, tolerance = 1e-3)
  
  problem <- clear_solution(problem)
  result <- solve(problem, SOLVER, qcp = TRUE, low = 12)
  expect_equal(result$value, 12.0, tolerance = 1e-3)
  expect_equal(result$getValue(x), 12.0, tolerance = 1e-3)
  
  problem <- clear_solution(problem)
  result <- solve(problem, SOLVER, qcp = TRUE, low = 0, high = 100)
  expect_equal(result$value, 12.0, tolerance = 1e-3)
  expect_equal(result$getValue(x), 12.0, tolerance = 1e-3)
})

test_that("test basic maximization with interval", {
  x <- Variable()
  expr <- ceil(x)
  
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconvex(expr))
  expect_true(is_quasiconcave(expr))
  expect_false(is_convex(expr))
  expect_false(is_concave(expr))
  expect_false(is_dcp(expr))
  expect_false(is_dgp(expr))
  
  problem <- Problem(Maximize(expr), list(x >= 12, x <= 17))
  expect_true(is_dqcp(problem))
  expect_false(is_dcp(problem))
  expect_false(is_dgp(problem))
  
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(x), 17.0, tolerance = 1e-3)
})

test_that("test basic maximum", {
  xy <- Variable(2)
  x <- xy[1]
  y <- xy[2]
  expr <- max_elemwise(ceil(x), ceil(y))
  
  problem <- Problem(Minimize(expr), list(x >= 12, x <= 17, y >= 17.4))
  expect_true(is_dqcp(problem))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), 18.0)
  expect_lt(result$getValue(x), 17.1)
  expect_gt(result$getValue(x), 11.9)
  expect_gt(result$getValue(y), 17.3)
})

test_that("test basic minimum", {
  xy <- Variable(2)
  x <- xy[1]
  y <- xy[2]
  expr <- min_elemwise(ceil(x), ceil(y))
  
  problem <- Problem(Minimize(expr), list(x >= 11.9, x <= 15.8, y >= 17.4))
  expect_true(is_dqcp(problem))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), 16.0)
  expect_lt(result$getValue(x), 16.0)
  expect_gt(result$getValue(x), 14.9)
  expect_gt(result$getValue(y), 17.3)
})

test_that("test basic composition", {
  xy <- Variable(2)
  x <- xy[1]
  y <- xy[2]
  expr <- max_elemwise(ceil(ceil(x)), ceil(ceil(y)))
  
  problem <- Problem(Minimize(expr), list(x >= 12, x <= 17, y >= 17.4))
  expect_true(is_dqcp(problem))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), 18.0)
  expect_lt(result$getValue(x), 17.1)
  expect_gt(result$getValue(x), 11.9)
  expect_gt(result$getValue(y), 17.3)
  
  # This problem should have the same solution.
  expr <- max_elemwise(floor(ceil(x)), floor(ceil(y)))
  problem <- Problem(Minimize(expr), list(x >= 12, x <= 17, y >= 17.4))
  expect_true(is_dqcp(problem))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), 18.0)
  expect_lt(result$getValue(x), 17.1)
  expect_gt(result$getValue(x), 11.9)
  result_gt(result$getValue(y), 17.3)
})

test_that("test basic floor", {
  x <- Variable()
  expr <- floor(x)
  
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconvex(expr))
  expect_true(is_quasiconcave(expr))
  expect_false(is_convex(expr))
  expect_false(is_concave(expr))
  expect_false(is_dcp(expr))
  expect_false(is_dgp(expr))
  
  problem <- Problem(Minimize(expr), list(x >= 11.8, x <= 17))
  expect_true(is_dqcp(problem))
  expect_false(is_dcp(problem))
  expect_false(is_dgp(problem))
  
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), 11.0)
  expect_gt(result$getValue(x), 11.7)
})

test_that("test basic multiply nonneg", {
  xy <- Variable(2, nonneg = TRUE)
  x <- xy[1]
  y <- xy[2]
  expr <- x*y
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconcave(expr))
  expect_false(is_quasiconvex(expr))
  
  expect_false(is_dcp(expr))
  
  problem <- Problem(Maximize(expr), list(x <= 12, y <= 6))
  expect_true(is_dqcp(problem))
  expect_false(is_dcp(problem))
  expect_false(is_dgp(problem))
  
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), 72, tolerance = 1e-1)
  expect_equal(result$getValue(x), 12, tolerance = 1e-1)
  expect_equal(result$getValue(y), 6, tolerance = 1e-1)
})

test_that("test basic multiply nonpos", {
  xy <- Variable(2, nonpos = TRUE)
  x <- xy[1]
  y <- xy[2]
  expr <- x*y
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconcave(expr))
  expect_false(is_quasiconvex(expr))
  
  expect_false(is_dcp(expr))
  
  problem <- Problem(Maximize(expr), list(x >= -12, y >= -6))
  expect_true(is_dqcp(problem))
  expect_false(is_dcp(problem))
  expect_false(is_dgp(problem))
  
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), 72, tolerance = 1e-1)
  expect_equal(result$getValue(x), -12, tolerance = 1e-1)
  expect_equal(result$getValue(y), -6, tolerance = 1e-1)
})

test_that("test basic multiply qcvx", {
  x <- Variable(nonneg = TRUE)
  y <- Variable(nonpos = TRUE)
  expr <- x*y
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconvex(expr))
  expect_false(is_quasiconcave(expr))
  
  expect_false(is_dcp(expr))
  
  problem <- Problem(Minimize(expr), list(x <= 7, y >= -6))
  expect_true(is_dqcp(problem))
  expect_false(is_dcp(problem))
  expect_false(is_dgp(problem))
  
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), -42, tolerance = 1e-1)
  expect_equal(result$getValue(x), 7, tolerance = 1e-1)
  expect_equal(result$getValue(y), -6, tolerance = 1e-1)
  
  x <- Variable(nonneg = TRUE)
  y <- Variable(nonpos = TRUE)
  expr <- y*x
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconvex(expr))
  expect_false(is_quasiconcave(expr))
  
  expect_false(is_dcp(expr))
  
  problem <- Problem(Minimize(expr), list(x <= 7, y >= -6))
  expect_true(is_dqcp(problem))
  expect_false(is_dcp(problem))
  expect_false(is_dgp(problem))
  
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), -42, tolerance = 1e-1)
  expect_equal(result$getValue(x), 7, tolerance = 1e-1)
  expect_equal(result$getValue(y), -6, tolerance = 1e-1)
})

test_that("test concave multiply", {
  xy <- Variable(2, nonneg = TRUE)
  x <- xy[1]
  y <- xy[2]
  expr <- sqrt(x)*sqrt(y)
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconcave(expr))
  expect_false(is_quasiconvex(expr))
  
  problem <- Problem(Maximize(expr), list(x <= 4, y <= 9))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), 6, tolerance = 1e-1)
  expect_equal(result$getValue(x), 4, tolerance = 1e-1)
  expect_equal(result$getValue(y), 9, tolerance = 1e-1)
  
  xy <- Variable(2, nonneg = TRUE)
  x <- xy[1]
  y <- xy[2]
  expr <- (sqrt(x) + 2.0)*(sqrt(y) + 4.0)
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconcave(expr))
  expect_false(is_quasiconvex(expr))
  
  problem <- Problem(Maximize(expr), list(x <= 4, y <= 9))
  result <- solve(problem, SOLVER, qcp = TRUE)
  # (2 + 2)*(3 + 4) = 28.
  expect_equal(result$getValue(problem@objective), 28, tolerance = 1e-1)
  expect_equal(result$getValue(x), 4, tolerance = 1e-1)
  expect_equal(result$getValue(y), 9, tolerance = 1e-1)
})

test_that("test basic ratio", {
  x <- Variable()
  y <- Variable(nonneg = TRUE)
  expr <- x/y
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconcave(expr))
  expect_true(is_quasiconvex(expr))
  
  problem <- Problem(Minimize(expr), list(x == 12, y <= 6))
  expect_true(is_dqcp(problem))
  
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), 2.0, tolerance = 1e-1)
  expect_equal(result$getValue(x), 12, tolerance = 1e-1)
  expect_equal(result$getValue(y), 6, tolerance = 1e-1)
  
  x <- Variable()
  y <- Variable(nonpos = TRUE)
  expr <- x/y
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconcave(expr))
  expect_true(is_quasiconvex(expr))
  
  problem <- Problem(Maximize(expr), list(x == 12, y >= -6))
  expect_true(is_dqcp(problem))
  
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), -2.0, tolerance = 1e-1)
  expect_equal(result$getValue(x), 12, tolerance = 1e-1)
  expect_equal(result$getValue(y), -6, tolerance = 1e-1)
})

test_that("test lin frac", {
  x <- Variable(2, nonneg = TRUE)
  A <- rbind(c(1.0, 2.0), c(3.0, 4.0))
  b <- c(0, 1)
  C <- 2*A
  d <- c(0, 1)
  lin_frac <- ((A %*% x) + b)/((C %*% x) + d)
  expect_true(is_dqcp(lin_frac))
  expect_true(is_quasiconvex(lin_frac))
  expect_true(is_quasiconcave(lin_frac))
  
  problem <- Problem(Minimize(sum(x)), list(x >= 0, lin_frac <= 1))
  expect_true(is_dqcp(problem))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), 0, tolerance = 1e-1)
  expect_equal(result$getValue(x), 0, tolerance = 1e-5)
})

test_that("test concave frac", {
  x <- Variable(nonneg = TRUE)
  concave_frac <- sqrt(x)/exp(x)
  expect_true(is_dqcp(concave_frac))
  expect_true(is_quasiconcave(concave_frac))
  expect_false(is_quasiconvex(concave_frac))
  
  problem <- Problem(Maximize(concave_frac))
  expect_true(is_dqcp(problem))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), 0.428, tolerance = 1e-1)
  expect_equal(result$getValue(x), 0.5, tolerance = 1e-1)
})

test_that("test length", {
  x <- Variable(5)
  expr <- length(x)
  expect_true(is_dqcp(expr))
  expect_true(is_quasiconvex(expr))
  expect_false(is_quasiconcave(expr))
  
  problem <- Problem(Minimize(expr), list(x[1] == 2.0, x[2] == 1.0))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), 2)
  expect_equal(result$getValue(x), matrix(c(2, 1, 0, 0, 0)), tolerance = 1e-7)
})

test_that("test length example", {
  n <- 10
  set.seed(1)
  A <- matrix(rnorm(n^2), nrow = n, ncol = n)
  x_star <- matrix(rnorm(n), nrow = n, ncol = 1)
  b <- A %*% x_star
  epsilon <- 1e-2
  x <- Variable(n)
  mse <- sum_squares(A %*% x - b)/n
  problem <- Problem(Minimize(length(x)), list(mse <= epsilon))
  if(!is_dqcp(problem))
    stop("Problem must be DQCP")
  # expect_true(is_dqcp(problem))
  
  result <- solve(problem, qcp = TRUE)
  if(is.na(result$value) || abs(result$value - 8) > (1e-8 + 1e-5*8))
    stop("Optimal value should be 8")
  # expect_equal(result$value, 8, tolerance = 1e-5)
})

test_that("test infeasible", {
  x <- Variable(2)
  problem <- Problem(Minimize(length(x)), list(x == -1, ceil(x) >= 1))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_in(result$status, c(INFEASIBLE, INFEASIBLE_INACCURATE))
})

test_that("test sign", {
  
})

test_that("test dist ratio", {
  
})

test_that("test infeasible exp constr", {
  
})

test_that("test infeasible inv_pos constr", {
  
})

test_that("test infeasible logistic constr", {
  
})

test_that("test noop exp constr", {
  
})

test_that("test noop inv_pos constr", {
  
})

test_that("test noop logistic constr", {
  
})

test_that("test gen_lambda_max matrix completion", {
  
})

test_that("test condition number", {
  
})

test_that("test card ls", {
  
})

test_that("test multiply const", {
  
})

test_that("test div const", {
  x <- Variable()
  obj <- Minimize(ceil(x)/0.5)
  problem <- Problem(obj, list(x >= 10))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(x), 10, tolerance = 1e-1)
  expect_equal(result$value, 20, tolerance = 1e-1)
  
  x <- Variable()
  obj <- Maximize(ceil(x)/-0.5)
})

test_that("test reciprocal", {
  x <- Variable(pos = TRUE)
  problem <- Problem(Minimize(1/x))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$value, 0, tolerance = 1e-3)
})

test_that("test abs", {
  x <- Variable(pos = TRUE)
  problem <- Problem(Minimize(abs(1/x)))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$value, 0, tolerance = 1e-3)
  
  x <- Variable(neg = TRUE)
  problem <- Problem(Minimize(abs(1/x)))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$value, 0, tolerance = 1e-3)
})

test_that("test tutorial example", {
  x <- Variable()
  y <- Variable(pos = TRUE)
  objective_fn <- -sqrt(x)/y
  problem <- Problem(Minimize(objective_fn), list(exp(x) <= y))
  # Smoke test.
  result <- solve(problem, SOLVER, qcp = TRUE)
})

test_that("test curvature", {
  x <- Variable(3)
  expr <- length(x)
  expect_equal(curvature(expr), QUASICONVEX)
  expr <- -length(x)
  expect_equal(curvature(expr), QUASICONCAVE)
  expr <- ceil(x)
  expect_equal(curvature(expr), QUASILINEAR)
  expect_true(is_quasilinear(expr))
})

test_that("test tutorial dqcp", {
  # The sign of variables affects curvature analysis.
  x <- Variable(nonneg = TRUE)
  concave_frac <- x*sqrt(x)
  constraint <- list(ceil(x) <= 10)
  problem <- Problem(Maximize(concave_frac), constraint)
  expect_true(is_quasiconcave(concave_frac))
  expect_true(is_dqcp(constraints[[1]]))
  expect_true(is_dqcp(problem))
  
  w <- Variable()
  fn <- w*sqrt(w)
  problem <- Problem(Maximize(fn))
  expect_false(is_dqcp(fn))
  expect_false(is_dqcp(problem))
})

test_that("test add constant", {
  # The sign of variables affects curvature analysis.
  x <- Variable()
  problem <- Problem(Minimize(ceil(x) + 5), list(x >= 2))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(x), 2, tolerance = 1e-7)
  expect_equal(result$getValue(problem@objective), 7, tolerance = 1e-7)
})

test_that("test max", {
  x <- Variable(2, pos = TRUE)
  obj <- max((1 - 2*sqrt(x) + x)/x)
  problem <- Problem(Minimize(obj), list(x[1] <= 0.5, x[2] <= 0.9))
  expect_true(is_dqcp(problem))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), 0.1715, tolerance = 1e-3)
})

test_that("test min", {
  x <- Variable(2)
  expr <- min(ceil(x))
  problem <- Problem(Maximize(expr), list(x[1] >= 11.9, x[1] <= 15.8, x[2] >= 17.4))
  expect_true(is_dqcp(problem))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$getValue(problem@objective), 16.0)
  expect_lt(result$getValue(x[1]), 16.0)
  expect_gt(result$getValue(x[1]), 14.9)
  expect_gt(result$getValue(x[2]), 17.3)
})

test_that("test sum of qccv not dqcp", {
  t <- Variable(5, pos = TRUE)
  expr <- sum(square(t)/t)
  expect_false(is_dqcp(expr))
})

test_that("test flip bounds", {
  x <- Variable(pos = TRUE)
  problem <- Problem(Maximize(ceil(x)), list(x <= 1))
  result <- solve(problem, SOLVER, qcp = TRUE, low = 0, high = 0.5)
  expect_gt(result$getValue(x), 0)
  expect_lte(result$getValue(x), 1)
  
  result <- solve(problem, SOLVER, qcp = TRUE, low = 0, high = NA)
  expect_gt(result$getValue(x), 0)
  expect_lte(result$getValue(x), 1)
  
  result <- solve(problem, SOLVER, qcp = TRUE, low = NA, high = 0.5)
  expect_gt(result$getValue(x), 0)
  expect_lte(result$getValue(x), 1)
})

test_that("test scalar sum", {
  x <- Variable(pos = TRUE)
  problem <- Problem(Minimize(sum(1/x)))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$value, 0, tolerance = 1e-3)
  
  problem <- Problem(Minimize(cumsum(1/x)))
  result <- solve(problem, SOLVER, qcp = TRUE)
  expect_equal(result$value, 0, tolerance = 1e-3)
})
