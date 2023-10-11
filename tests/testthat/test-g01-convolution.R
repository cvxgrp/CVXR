context("test-g01-convolution")
TOL <- 1e-6

test_that("test 1D convolution", {
  skip_on_cran()
  n <- 3
  x <- Variable(n)
  f <- c(1, 2, 3)
  g <- c(0, 1, 0.5)
  f_conv_g <- c(0, 1, 2.5, 4, 1.5)
  expr <- conv(f, g)
  expect_true(is_constant(expr))
  expect_equal(dim(expr), c(5, 1))
  expect_equal(value(expr), f_conv_g)

  expr <- conv(f, x)
  expect_true(is_affine(expr))
  expect_equal(dim(expr), c(5, 1))

  # Matrix stuffing
  prob <- Problem(Minimize(p_norm(expr, 1)), list(x == g))
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, sum(f_conv_g), tolerance = 1e-3)
  expect_equal(result$getValue(expr), f_conv_g, tolerance = TOL)
  
  # # Expression trees.
  # prob <- Problem(Minimize(p_norm(expr, 1)))
  # prob_mat_vs_mul_funcs(prob)
  # result <- solve(prob, solver = "SCS", expr_tree = TRUE, verbose =  TRUE)
  # expect_equal(result$value, 0, tolerance = 1e-1)
})

prob_mat_vs_mul_funcs <- fucntion(prob) {
  tmp <- get_problem_data(prob, solver = "SCS")
  data <- tmp[[1]]
  dims <- tmp[[2]]
  A <- data$A
  canon <- canonicalize(prob, "SCS")
  objective <- canon[[1]]
  constr_map <- canon[[2]]
  dims <- canon[[3]]
  solver <- canon[[4]]
  
  all_ineq <- c(constr_map[[EQ]], constr_map[[LEQ]])
  tmp <- .get_var_offsets(prob, objective, all_ineq)
  var_offsets <- tmp[[1]]
  var_sizes <- tmp[[2]]
  x_length <- tmp[[3]]
  
  constraints <- c(constr_map[[EQ]], constr_map[[LEQ]])
  constraints <- prune_constants(constraints)
  tmp <- get_mul_funcs(constraints, dims, var_offsets, var_sizes, x_length)
  Amul <- tmp[[1]]
  ATmul <- tmp[[2]]
  
  vec <- matrix(c(1:x_length))
  
  # A %*% vec.
  result <- rep(0, nrow(A))
  Amul(vec, result)
  expect_equal(A %*% vec, result, tolerance = TOL)
  Amul(vec, result)
  expect_equal(2*A %*% vec, result, tolerance = TOL)
  
  # t(A) %*% vec.
  vec <- matrix(0:(nrow(A) - 1))
  result <- rep(0, ncol(A))
  ATmul(vec, result)
  expect_equal(t(A) %*% vec, result, tolerance = TOL)
  ATmul(vec, result)
  expect_equal(2*t(A) %*% vec, result, tolerance = TOL)
}

mat_from_func <- function(func, rows, cols) {
  # Convert a multiplier function to a matrix.
  test_vec <- rep(0, cols)
  result <- rep(0, rows)
  mat <- matrix(0, nrow = rows, ncol = cols)
  for(i in 1:cols) {
    test_vec[i] <- 1.0
    func(test_vec, result)
    mat[,i] <- result
    test_vec <- test_vec*0
    result <- result*0
  }
  return(mat)
}

test_that("test a problem with convolution", {
  skip_on_cran()
  N <- 5
  y <- matrix(stats::rnorm(N), nrow = N, ncol = 1)
  h <- matrix(stats::rnorm(2), nrow = 2, ncol = 1)
  x <- Variable(N)
  v <- conv(h, x)
  obj <- Minimize(sum(multiply(y, v[1:N])))
  result <- solve(Problem(obj, list()), solver = "ECOS")
  print(result$value)
})

test_that("test convolution with 0D input", {
  x <- Variable(1)
  problem <- Problem(Minimize(max(conv(1, multiply(1, x)))), list(x >= 0))
  result <- solve(problem, solver = "ECOS")
  expect_equal(result$status, "optimal")
})
