context("test_gurobi_write")

test_that("test the Gurobi model write", {
  skip_on_cran()
  if(!("GUROBI" %in% installed_solvers())) {
    print("Gurobi is not installed. Skipping test.")
    return()
  }
  
  if(!dir.exists("./resources/"))
    dir.create("./resources/")
  else if(file.exists("./resources/gurobi_model.lp"))
    file.remove("./resources/gurobi_model.lp")
  
  m <- 20
  n <- 15
  set.seed(0)
  
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  b <- matrix(rnorm(m))
  
  x <- Variable(n)
  cost <- sum_squares(A %*% x - b)
  prob <- Problem(Minimize(cost))
  result <- solve(prob, solver = "GUROBI", verbose = TRUE, save_file = "./resources/gurobi_model.lp")
  expect_true(file.exists("./resources/gurobi_model.lp"))
})