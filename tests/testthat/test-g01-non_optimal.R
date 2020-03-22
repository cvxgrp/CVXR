context("test-g01-non_optimal")

test_that("Test scalar LP problems", {
  skip_on_cran()
  x1 <- Variable()
  x2 <- Variable()
  obj <- Minimize(-x1-x2)
  constraints <- list(2*x1 + x2 >= 1, x1 + 3*x2 >= 1, x1 >= 0, x2 >= 0)
  p_unb <- Problem(obj, constraints)
  p_inf <- Problem(Minimize(x1), list(x1 >= 0, x1 <= -1))
  for(solver in c("ECOS", "SCS")) {
    print(solver)
    result <- solve(p_unb, solver = solver)
    expect_equal(tolower(result$status), "unbounded")
    result <- solve(p_inf, solver = solver)
    expect_equal(tolower(result$status), "infeasible")
  }
})

test_that("Test vector LP problems", {
  skip_on_cran()
  # Infeasible and unbounded problems
  x <- Variable(5)
  p_inf <- Problem(Minimize(sum(x)), list(x >= 1, x <= 0))
  p_unb <- Problem(Minimize(sum(x)), list(x <= 1))
  for(solver in c("ECOS", "SCS")) {
    print(solver)
    result <- solve(p_unb, solver = solver)
    expect_equal(tolower(result$status), "unbounded")
    result <- solve(p_inf, solver = solver)
    expect_equal(tolower(result$status), "infeasible")
  }
})

# test_that("Test the optimal inaccurate status", {
#  x <- Variable(5)
#  prob <- Problem(Maximize(sum(sqrt(x))), list(x <= 0))
#  result <- solve(prob, solver = "SCS")
#  expect_equal(tolower(result$status), "optimal_inaccurate")
#  expect_false(is.na(result$value))
# })

# test_that("Test SOC problems", {
#   # Infeasible and unbounded problems.
#   x <- Variable(5)
#   obj <- Maximize(sum(sqrt(x)))
#   p_inf <- Problem(obj, list(x >= 1, x <= 0))
#   p_unb <- Problem(obj, list(x >= 1))
#   for(solver in c("ECOS", "SCS")) {
#     print(solver)
#     result <- solve(p_inf, solver = solver)
#     expect_equal(result$status, "infeasible")
#     result <- solve(p_unb, solver = solver)
#     expect_equal(result$status, "unbounded")
#   }
# })

# test_that("Test PSD problems", {
#   # Infeasible and unbounded problems.
#   X <- Variable(5,5)
#   obj <- Maximize(lambda_min(X))
#   p_inf <- Problem(obj, list(X >= 1, X <= 0))
#   p_unb <- Problem(obj)
#   for(solver in c("SCS")) {
#     print(solver)
#     result <- solve(p_inf, solver = solver)
#     expect_equal(result$status, "infeasible")
#     result <- solve(p_unb, solver = solver)
#     expect_equal(result$status, "unbounded")
#   }
# })
