TOL <- 1e-6

x_bool <- Bool()
y_int <- Int()
A_bool <- Bool(3,2)
B_int <- Int(2,3)

test_that("Test that MIP problems are deterministic", {
  data_recs <- list()
  result_recs <- list()
  for(i in 1:5) {
    obj <- Minimize(square(y_int - 0.2))
    p <- Problem(obj, list(A_bool == 0, x_bool == B_int))
    data_recs <- c(data_recs, list(get_problem_data(p, "ECOS_BB")))
  }
  
  # Check that problem data and result is always the same
  for(i in 1:5) {
    for(key in c("c", "A", "b", "G", "h", "bool_vars_idx", "int_vars_idx")) {
      lh_item <- data_recs[[1]][[key]]
      rh_item <- data_recs[[i]][[key]]
      if(key %in% c("A", "G")) {
        lh_item <- as.matrix(lh_item)
        rh_item <- as.matrix(rh_item)
      }
      expect_equal(lh_item, rh_item, tolerance = TOL)
    }
  }
})

test_that("Test Boolean problems", {
  # Bool in objective
  obj <- Minimize(square(x_bool - 0.2))
  p <- Problem(obj, list())
  result <- solve(p)
  expect_equal(result$value, 0.04, tolerance = TOL)
  expect_equal(result$getValue(x_bool), 0, tolerance = TOL)
  
  # Bool in constraint
  t <- Variable()
  obj <- Minimize(t)
  p <- Problem(obj, list(square(x_bool) <= t))
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

test_that("Test Integer problems", {
  # Int in objective
  obj <- Minimize(square(y_int - 0.2))
  p <- Problem(obj, list())
  result <- solve(p)
  expect_equal(result$value, 0.04, tolerance = TOL)
  expect_equal(result$getValue(y_int), 0, tolerance = TOL)
})
