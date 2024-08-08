context("ram_limited")
TOL <- 1e-^6

test_that("RAM Limited: issue826", {
  # In cvxpy GitHub issue #826, it was discovered that cvxcore's C++
  # implementation implicitly limited problem data (such as a
  # constraint matrix) to have at most 2^(32)-1 nonzero entries.
  #
  # This test is for checking that #826 is resolved.
  n <- 2^8
  m <- as.integer(2^32/n) + 1
  
  vals <- seq(0, m*n - 1)/1000
  A <- matrix(vals, nrow = n, ncol = m, byrow = TRUE)
  x <- Variable(m)
  cons <- list(A %*% x >= 0)
  prob <- Problem(Maximize(0), cons)
  data <- get_problem_data(prob, solver = "SCS")
  vals_canon <- data[[1]]$A@data
  
  vdiff <- vals - vals_canon
  verr <- abs(vdiff)
  expect_lte(verr, 1e-3)
  print("issue826 test finished")
})