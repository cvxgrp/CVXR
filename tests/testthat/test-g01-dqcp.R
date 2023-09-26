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
  
})

test_that("test concave multiply", {
  
})

test_that("test basic ratio", {
  
})

test_that("test lin frac", {
  
})

test_that("test concave frac", {
  
})

test_that("test length", {
  
})

test_that("test length example", {
  
})

test_that("test infeasible", {
  
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
  
})

test_that("test reciprocal", {
  
})

test_that("test abs", {
  
})

test_that("test tutorial example", {
  
})

test_that("test curvature", {
  
})

test_that("test tutorial dqcp", {
  
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
