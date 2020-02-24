## These tests are problematic on CRAN because they
## actually install RMOSEK and Rcplex on debian and
## they never tell us what IBM cplex is being used.
## So there is no real way to troubleshoot, short of
## installing everything ourselves.

test_that("Test Boolean problems", {
  skip_on_cran()
  # Bool in objective
  obj <- Minimize((x_bool - 0.2)^2)
  p <- Problem(obj, list())
  result <- solve(p)
  expect_equal(result$value, 0.04, tolerance = TOL)
  expect_equal(result$getValue(x_bool), 0, tolerance = TOL)

  # Bool in constraint
  t <- Variable()
  obj <- Minimize(t)
  p <- Problem(obj, list(x_bool^2 <= t))
  result <- solve(p)
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(x_bool), 0, tolerance = 1e-4)

  # Matrix Bool in objective
  C <- cbind(c(0,1,0), c(1,1,1))
  obj <- Minimize(sum_squares(A_bool - C))
  p <- Problem(obj, list())
  result <- solve(p)
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(A_bool), C, tolerance = 1e-4)

  # Matrix Bool in constraint
  t <- Variable()
  obj <- Minimize(t)
  p <- Problem(obj, list(sum_squares(A_bool - C) <= t))
  result <- solve(p)
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(A_bool), C, tolerance = 1e-4)
})
