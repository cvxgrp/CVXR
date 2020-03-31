context("test-g01-mip_vars")
TOL <- 1e-6

x_bool <- Variable(boolean = TRUE)
y_int <- Variable(integer = TRUE)
A_bool <- Variable(3, 2, boolean = TRUE)
B_int <- Variable(2, 3, integer = TRUE)

test_that("Test that MIP problems are deterministic", {
  skip_on_cran()
  data_recs <- list()
  result_recs <- list()
  for(i in 1:5) {
    obj <- Minimize((y_int - 0.2)^2)
    p <- Problem(obj, list(A_bool == 0, x_bool == B_int))
    data_recs <- c(data_recs, list(get_problem_data(p, "ECOS_BB")))
  }

  # Check that problem data and result is always the same
  for(i in 1:5) {
    for(key in c("c", "A", "b", "G", "h", "bool_vars_idx", "int_vars_idx")) {
      lh_item <- data_recs[[1]][[1]][[key]]
      rh_item <- data_recs[[i]][[1]][[key]]
      if(key %in% c("A", "G")) {
        lh_item <- as.matrix(lh_item)
        rh_item <- as.matrix(rh_item)
      }
      expect_equal(lh_item, rh_item, tolerance = TOL)
    }
  }
})

test_that("Test Boolean problems", {
  skip_on_cran()
  # Bool in objective
  obj <- Minimize((x_bool - 0.2)^2)
  p <- Problem(obj, list())
  result <- solve(p, solver = "ECOS_BB")
  expect_equal(result$value, 0.04, tolerance = TOL)
  expect_equal(result$getValue(x_bool), 0, tolerance = TOL)

  # Bool in constraint
  t <- Variable()
  obj <- Minimize(t)
  p <- Problem(obj, list(x_bool^2 <= t))
  result <- solve(p, solver = "ECOS_BB")
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(x_bool), 0, tolerance = 1e-4)

  # Matrix Bool in objective
  C <- cbind(c(0,1,0), c(1,1,1))
  obj <- Minimize(sum_squares(A_bool - C))
  p <- Problem(obj, list())
  result <- solve(p, solver = "ECOS_BB")
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(A_bool), C, tolerance = 1e-4)

  # Matrix Bool in constraint
  t <- Variable()
  obj <- Minimize(t)
  p <- Problem(obj, list(sum_squares(A_bool - C) <= t))
  result <- solve(p, solver = "ECOS_BB")
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(A_bool), C, tolerance = 1e-4)
})

test_that("Test Integer problems", {
  skip_on_cran()
  # Int in objective
  obj <- Minimize((y_int - 0.2)^2)
  p <- Problem(obj, list())
  result <- solve(p, solver = "ECOS_BB")
  expect_equal(result$value, 0.04, tolerance = TOL)
  expect_equal(result$getValue(y_int), 0, tolerance = TOL)

  # Infeasible integer problem
  obj <- Minimize(0)
  p <- Problem(obj, list(y_int == 0.5))
  result <- solve(p, solver = "ECOS_BB")   # TODO_NARAS_10: Returns solver_error (MOSEK) or freezes (ECOS_BB). DONE by Naras
  ##result <- solve(p, solver = "CBC")
  expect_true(grepl("infeasible", result$status))
})

test_that("Test SOCP problems", {
  skip_on_cran()
  # Int in objective
  t <- Variable()
  obj <- Minimize(t)
  p <- Problem(obj, list((y_int - 0.2)^2 <= t))
  result <- solve(p, solver = "ECOS_BB")
  expect_equal(result$value, 0.04, tolerance = TOL)
  expect_equal(result$getValue(y_int), 0, tolerance = 2 * TOL)
})

test_that("Test Boolean SOCP problems", {
  skip_on_cran()
  # Int in objective
  t <- Variable()
  obj <- Minimize(t)
  p <- Problem(obj, list((x_bool - 0.2)^2 <= t))
  result <- solve(p, solver = "ECOS_BB")
  expect_equal(result$value, 0.04, tolerance = TOL)
  expect_equal(result$getValue(x_bool), 0, tolerance = TOL)
})
