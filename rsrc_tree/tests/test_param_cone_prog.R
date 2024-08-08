context("test_param_cone_prog")
TOL <- 1e-2

test_that("test log problem", {
  # Log in objective.
  x <- Variable(2)
  x_id_char <- as.character(id(x))
  var_dict <- list()
  var_dict[[x_id_char]] <- x
  obj <- Maximize(sum(log(x)))
  constr <- list(x <= c(1, exp(1)))
  problem <- Problem(obj, constr)
  data <- get_problem_data(problem, solver = "SCS")[[1]]
  param_cone_prog <- data[[PARAM_PROB]]
  solver <- SCS()
  raw_solution <- solve_via_data(solver, data, warm_start = FALSE, verbose = FALSE, solver_opts = list())$x
  sltn_dict <- split_solution(param_cone_prog, raw_solution, active_vars = var_dict)
  adjoint <- split_adjoint(param_cone_prog, sltn_dict)
  expect_equal(dim(adjoint), dim(raw_solution))
  for(value in sltn_dict[[x_id_char]])
    expect_true(any(value == adjoint))
  
  result <- solve(problem, solver = "SCS")
  expect_equal(result$getValue(x), sltn_dict[[x_id_char]], tolerance = TOL)
  
  # Log in constraint.
  obj <- Minimize(sum(x))
  constr <- list(log(x) >= 0, x <= c(1,1))
  problem <- Problem(obj, constr)
  data <- get_problem_data(problem, solver = "SCS")[[1]]
  param_cone_prog <- data[[PARAM_PROB]]
  solver <- SCS()
  raw_solution <- solve_via_data(solver, data, warm_start = FALSE, verbose = FALSE, solver_opts = list())$x
  sltn_dict <- split_solution(param_cone_prog, raw_solution, active_vars = var_dict)
  adjoint <- split_adjoint(param_cone_prog, sltn_dict)
  expect_equal(dim(adjoint), dim(raw_solution))
  for(value in sltn_dict[[x_id_char]])
    expect_true(any(value == adjoint))
  expect_equal(split_adjoint(param_cone_prog, sltn_dict), raw_solution, tolerance = TOL)
  
  result <- solve(problem, solver = "SCS")
  expect_equal(result$getValue(x), sltn_dict[[x_id_char]], tolerance = TOL)
})

test_that("test PSD variable", {
  s <- Variable(2, 2, PSD = TRUE)
  s_id_char <- as.character(id(s))
  var_dict <- list()
  var_dict[[s_id_char]] <- s
  obj <- Maximize(min_elemwise(s[1,2], 10))
  const <- list(diag(s) == c(1,1))
  problem <- Problem(obj, const)
  data <- get_problem_data(problem, solver = "SCS")[[1]]
  param_cone_prog <- data[[PARAM_PROB]]
  solver <- SCS()
  raw_solution <- solve_via_data(solver, data, warm_start = FALSE, verbose = FALSE, solver_opts = list())$x
  sltn_dict <- split_solution(param_cone_prog, raw_solution, active_vars = var_dict)
  expect_equal(dim(sltn_dict[[s_id_char]]), dim(s))
  sltn_value <- sltn_dict[[s_id_char]]
  adjoint <- split_adjoint(param_cone_prog, sltn_dict)
  expect_equal(dim(adjoint), dim(raw_solution))
  expect_true(any(sltn_value[1,1] == adjoint))
  expect_true(any(sltn_value[2,2] == adjoint))
  
  # Off-diagonals of adjoint will be scaled by two.
  expect_true(any(np.isclose(2*sltn_value[1,2], adjoint)))
  expect_true(any(np.isclose(2*sltn_value[2,1], adjoint)))
  
  result <- solve(problem, solver = "SCS")
  expect_equal(result$getValue(s), sltn_value, tolerance = TOL)
})
