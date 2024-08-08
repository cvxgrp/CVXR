context("test_valinvec2mixedint")

make_test_1 <- function(ineq_form) {
  # vec contains a contiguous range of integers.
  x <- Variable(4)
  expect_x <- matrix(c(0, 7, 3, 0))
  vec <- 0:9
  objective <- Maximize(x[1] + x[2] + 2*x[3] - 2*x[4])
  
  constr1 <- FiniteSet(x[1], vec, ineq_form = ineq_form)
  constr2 <- FiniteSet(x[2], vec, ineq_form = ineq_form)
  constr3 <- FiniteSet(x[3], vec, ineq_form = ineq_form)
  constr4 <- FiniteSet(x[4], vec, ineq_form = ineq_form)
  constr5 <- x[1] + 2*x[3] <= 700
  constr6 <- 2*x[2] - 8*x[3] <= 0
  constr7 <- x[2] - 2*x[3] + x[4] >= 1
  constr8 <- x[1] + x[2] + x[3] + x[4] == 10
  
  obj_pair <- list(objective, 13.0)
  con_pairs <- list(list(constr1, NA_real_),
                    list(constr2, NA_real_),
                    list(constr3, NA_real_),
                    list(constr4, NA_real_),
                    list(constr5, NA_real_),
                    list(constr6, NA_real_),
                    list(constr7, NA_real_),
                    list(constr8, NA_real_))
  var_pairs <- list(list(x, expect_x))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

make_test_2 <- function(ineq_form) {
  x <- Variable()
  expect_x <- matrix(-1.125)
  objective <- Minimize(x)
  vec <- c(-1.125, 1, 2)
  constr1 <- x >= -1.25
  constr2 <- x <= 10
  constr3 <- FiniteSet(x, vec, ineq_form = ineq_form)
  obj_pairs <- list(objective, -1.125)
  var_pairs <- list(list(x, expect_x))
  con_pairs <- list(list(constr1, NA_real_),
                    list(constr2, NA_real_),
                    list(constr3, NA_real_))
  sth <- SolverTestHelper(obj_pairs, var_pairs, con_pairs)
  return(sth)
}

make_test_3 <- function(ineq_form) {
  # Case when size(vec) == 1.
  x <- Variable()
  objective <- Minimize(abs(x-3))
  vec <- c(1)
  cons1 <- FiniteSet(x, vec, ineq_form = ineq_form)
  expected_x <- matrix(1)
  obj_pair <- list(objective, 2.0)
  var_pairs <- list(list(x, expected_x))
  con_pairs <- list(list(cons1, NA_real_))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

make_test_4 <- function(ineq_form) {
  # Case when vec houses duplicates.
  x <- Variable()
  objective <- Minimize(abs(x-3))
  vec <- c(1,1,1,2,2,3,3)
  cons1 <- FiniteSet(x, vec, ineq_form = ineq_form)
  expected_x <- matrix(3)
  obj_pair <- list(objective, 0.0)
  var_pairs <- list(list(x, expected_x))
  con_pairs <- list(list(cons1, NA_real_))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

make_test_5 <- function(ineq_form) {
  # Case when input expression to FiniteSet constraint is affine.
  x <- Variable(4)
  vec <- 0:9
  objective <- Maximize(x[1] + x[2] + 2*x[3] - 2*x[4])
  expr0 <- 2*x[1] + 1
  expr2 <- 3*x[3] + 5
  
  constr1 <- FiniteSet(expr0, vec, ineq_form = ineq_form)
  constr2 <- FiniteSet(x[2], vec, ineq_form = ineq_form)
  constr3 <- FiniteSet(expr2, vec, ineq_form = ineq_form)
  constr4 <- FiniteSet(x[4], vec, ineq_form = ineq_form)
  constr5 <- x[1] + 2*x[3] <= 700
  constr6 <- 2*x[2] - 8*x[3] <= 0
  constr7 <- x[2] - 2*x[3] + x[4] >= 1
  constr8 <- x[1] + x[2] + x[3] + x[4] == 10
  
  expected_x <- matrix(c(4,4,1,1))
  obj_pair <- list(objective, 8.0)
  con_pairs <- list(list(constr1, NA_real_),
                    list(constr2, NA_real_),
                    list(constr3, NA_real_),
                    list(constr4, NA_real_),
                    list(constr5, NA_real_),
                    list(constr6, NA_real_),
                    list(constr7, NA_real_),
                    list(constr8, NA_real_))
  var_pairs <- list(list(x, expected_x))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

make_test_6 <- function(ineq_form) {
  # vec contains only real quantities + passed expression is affine.
  x <- Variable()
  expect_x <- matrix(-1.0625)
  objective <- Minimize(x)
  vec <- c(-1.125, 1.5, 2.24)
  constr1 <- x >= -1.25
  constr2 <- x <= 10
  expr <- 2*x + 1
  constr3 <- FiniteSet(expr, vec, ineq_form = ineq_form)
  obj_pairs <- list(objective, -1.0625)
  var_pairs <- list(list(x, expect_x))
  con_pairs <- list(list(constr1, NA_real_),
                    list(constr2, NA_real_),
                    list(constr3, NA_real_))
  sth <- SolverTestHelper(obj_pairs, var_pairs, con_pairs)
  return(sth)
}

make_test_7 <- function(ineq_form) {
  # For testing vectorization of FiniteSet class.
  x <- Variable(4)
  expect_x <- matrix(c(0, 7, 3, 0))
  vec <- 0:9
  objective <- Maximize(x[1] + x[2] + 2*x[3] - 2*x[4])
  
  constr1 <- FiniteSet(x, vec, ineq_form = ineq_form)
  constr2 <- x[1] + 2*x[3] <= 700
  constr3 <- 2*x[2] - 8*x[3] <= 0
  constr4 <- x[2] - 2*x[3] + x[4] >= 1
  constr5 <- x[1] + x[2] + x[3] + x[4] == 10
  
  obj_pair <- list(objective, 13.0)
  con_pairs <- list(list(constr1, NA_real_),
                    list(constr2, NA_real_),
                    list(constr3, NA_real_),
                    list(constr4, NA_real_),
                    list(constr5, NA_real_))
  var_pairs <- list(list(x, expect_x))
  sth <- SolverTestHelper(obj_pair, var_pairs, con_pairs)
  return(sth)
}

TestFiniteSet.test_1 <- function(ineq_form) {
  sth <- make_test_1(ineq_form)
  result <- solve(sth, solver = "GLPK_MI")
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
}

TestFiniteSet.test_2 <- function(ineq_form) {
  sth <- make_test_2(ineq_form)
  result <- solve(sth, solver = "GLPK_MI")
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
}

TestFiniteSet.test_3 <- function(ineq_form) {
  sth <- make_test_3(ineq_form)
  result <- solve(sth, solver = "GLPK_MI")
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
}

TestFiniteSet.test_4 <- function(ineq_form) {
  sth <- make_test_4(ineq_form)
  result <- solve(sth, solver = "GLPK_MI")
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
}

TestFiniteSet.test_5 <- function(ineq_form) {
  sth <- make_test_5(ineq_form)
  result <- solve(sth, solver = "GLPK_MI")
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
}

TestFiniteSet.test_6 <- function(ineq_form) {
  sth <- make_test_6(ineq_form)
  result <- solve(sth, solver = "GLPK_MI")
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
}

TestFiniteSet.test_7 <- function(ineq_form) {
  sth <- make_test_7(ineq_form)
  result <- solve(sth, solver = "GLPK_MI")
  verify_objective(sth, result, tolerance = 1e-3)
  verify_primal_values(sth, result, tolerance = 1e-3)
}

TestFiniteSet.test_8 <- function(ineq_form) {
  # Test parametrized FiniteSet.
  x <- Variable()
  objective <- Maximize(x)
  set_vals <- Parameter(5, value = 0:4)
  constraints <- list(FiniteSet(x, set_vals, ineq_form = ineq_form))
  problem <- Problem(objective, constraints)
  result <- solve(problem, solver = "GLPK_MI")
  expect_true(is.allclose(result$getValue(x), 4))
}

TestFiniteSet.test_9 <- function(ineq_form) {
  # Test passing a set data structure from the R sets library.
  if(!require(sets)) {
    print("No R 'sets' library detected. Skipping test.")
    return()
  }
  
  x <- Variable()
  objective <- Maximize(x)
  set_vals <- sets::set(0:4)
  constraints <- list(FiniteSet(x, set_vals, ineq_form = ineq_form))
  problem <- Problem(objective, constraints)
  result <- solve(problem, solver = "GLPK_MI")
  expect_true(is.allclose(result$getValue(x), 4))
}

TestFiniteSet.test_10 <- function(ineq_form) {
  # Test set with two elements.
  x <- Variable()
  objective <- Maximize(x)
  set_vals <- list(1, 2)
  constraints <- list(FiniteSet(x, set_vals, ineq_form = ineq_form))
  problem <- Problem(objective, constraints)
  result <- solve(problem, solver = "GLPK_MI")
  expect_true(is.allclose(result$getValue(x), 2))
}

TestFiniteSet.test_11 <- function(ineq_form) {
  # Test 2D Variable.
  x <- Variable(2, 2)
  objective <- Maximize(sum(x))
  set_vals <- list(1, 2, 3)
  constraints <- list(FiniteSet(x, set_vals, ineq_form = ineq_form))
  problem <- Problem(objective, constraints)
  result <- solve(problem, solver = "GLPK_MI")
  expect_true(is.allclose(result$getValue(x), matrix(1, nrow = 2, ncol = 2)*max(unlist(set_vals))))
}

TestFiniteSet.test_non_affine_exception <- function(ineq_form) {
  # Exception test: non-affine expression.
  x <- Variable()
  x_abs <- abs(x)
  set_vals <- list(1, 2, 3)
  expect_error(FiniteSet(x_abs, set_vals, ineq_form = ineq_form), "must be affine")
}

TestFiniteSet.test_independent_entries <- function(ineq_form) {
  x <- Variable(2, 2)
  objective <- Maximize(sum(x))
  set_vals <- list(0, 1, 2)
  constraints <- list(FiniteSet(x, set_vals, ineq_form = ineq_form),
                      x <= matrix(0:4, nrow = 2, ncol = 2, byrow = TRUE))
  problem <- Problem(objective, constraints)
  result <- solve(problem, solver = "GLPK_MI")
  expect_true(is.allclose(result$getValue(x), rbind(c(0,1), c(2,2))))
}

test_that("test default argument", {
  if(!require(sets)) {
    print("No R 'sets' library detected. Skipping test.")
    return()
  }
  
  if(!("GLPK_MI" %in% installed_solvers())) {
    print("GLPK_MI solver is required, but not installed. Skipping test.")
    return()
  }
  
  x <- Variable()
  objective <- Maximize(x)
  set_vals <- sets::set(0:4)
  constraints <- list(FiniteSet(x, set_vals))
  problem <- Problem(objective, constraints)
  result <- solve(problem, solver = "GLPK_MI")
  expect_true(is.allclose(result$getValue(x), 4))
})

test_that("test all finite sets", {
  if(!("GLPK_MI" %in% installed_solvers())) {
    print("GLPK_MI solver is required, but not installed. Skipping test.")
    return()
  }
  
  for(ineq_form in c(TRUE, FALSE)) {
    TestFiniteSet.test_1(ineq_form)
    TestFiniteSet.test_2(ineq_form)
    TestFiniteSet.test_3(ineq_form)
    TestFiniteSet.test_4(ineq_form)
    TestFiniteSet.test_5(ineq_form)
    TestFiniteSet.test_6(ineq_form)
    TestFiniteSet.test_7(ineq_form)
    TestFiniteSet.test_8(ineq_form)
    TestFiniteSet.test_9(ineq_form)
    TestFiniteSet.test_10(ineq_form)
    TestFiniteSet.test_11(ineq_form)
    TestFiniteSet.test_non_affine_exception(ineq_form)
    TestFiniteSet.test_independent_entries(ineq_form)
  }
})
