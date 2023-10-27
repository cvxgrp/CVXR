context("test-g01-problem")
TOL <- 1e-6

a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

ConicSolver <- CVXR:::ConicSolver
INSTALLED_SOLVERS <- installed_solvers()
SCS.dims_to_solver_dict <- CVXR:::SCS.dims_to_solver_dict
ECOS.dims_to_solver_dict <- CVXR:::ECOS.dims_to_solver_dict
.p_norm <- CVXR:::.p_norm

test_that("test string representation", {
  skip_on_cran()
  obj <- Minimize(a)
  prob <- Problem(obj)
  expect_equal(as.character(prob), paste("Problem(", as.character(obj), ", ())", sep = ""))
  
  constraints <- list(x*2 == x, x == 0)
  prob <- Problem(obj, constraints)
  expect_equal(as.character(prob), "Problem(", as.character(obj), ", (", paste(sapply(constraints, as.character), sep = ", "), "))", sep = "")
  
  # Test str.
  a_name <- name(a)
  result <- paste("minimize ", a_name, "\nsubject to ", a_name, " == 0\n           ", name(a), " <= 0", sep = "")
  prob <- Problem(Minimize(a), list(ZeroConstraint(a), NonPosConstraint(a)))
  expect_equal(as.character(prob), result)
})

test_that("test the variables method", {
  skip_on_cran()
  p <- Problem(Minimize(a), list(a <= x, b <= A + 2))
  vars_ <- variables(p)
  ref <- list(a, x, b, A)
  expect_setequal(vars_, ref)
})

test_that("test var dict", {
  skip_on_cran()
  p <- Problem(Minimize(a), list(a <= x, b <= A + 2))
  expect_equal(p@var_dict, list(a = a, x = x, b = b, A = A))
})

test_that("test the parameters method", {
  skip_on_cran()
  p1 <- Parameter()
  p2 <- Parameter(3, nonpos = TRUE)
  p3 <- Parameter(4, 4, nonneg = TRUE)
  p <- Problem(Minimize(p1), list(a + p1 <= p2, b <= p3 + p3 + 2))
  params <- parameters(p)
  ref <- c(p1, p2, p3)
  expect_setequal(params, ref)
})

test_that("test param dict", {
  skip_on_cran()
  p1 <- Parameter(name = "p1")
  p2 <- Parameter(3, nonpos = TRUE, name = "p2")
  p3 <- Parameter(4, 4, nonneg = TRUE, name = "p3")
  p <- Problem(Minimize(p1), list(a + p1 <= p2, b <= p3 + p3 + 2))
  expect_equal(p@param_dict, list(p1 = p1, p2 = p2, p3 = p3))
})

test_that("test solving a problem with unspecified parameters", {
  skip_on_cran()
  param <- Parameter(name = "lambda")
  problem <- Problem(Minimize(param), list())
  expect_error(solve(problem, solver = "SCS"), "A Parameter (whose name is 'lambda')")
})

test_that("test the constants method", {
  skip_on_cran()
  c1 <- matrix(stats::rnorm(2), nrow = 1, ncol = 2)
  c2 <- matrix(stats::rnorm(2), nrow = 2, ncol = 1)
  p <- Problem(Minimize(c1 %*% x), list(x >= c2))
  constants_ <- constants(p)
  ref <- list(c1, c2)
  expect_equal(length(constants_), length(ref))
  mapply(function(c, r) {
    expect_equal(dim(c), dim(r))
    expect_true(all(value(c) == r))
  }, constants_, ref)

  # Single scalar constants.
  p <- Problem(Minimize(a), list(x >= 1))
  constants_ <- constants(p)
  ref <- list(matrix(1))
  expect_equal(length(ref), length(constants_))
  mapply(function(c, r) {
    expect_equal(dim(c), dim(r))
    expect_true(all(value(c) == r))
  }, constants_, ref)
})

test_that("Test the size_metrics method", {
  skip_on_cran()
  p1 <- Parameter()
  p2 <- Parameter(3, nonpos = TRUE)
  p3 <- Parameter(4, 4, nonneg = TRUE)

  c1 <- matrix(stats::rnorm(2), nrow = 2, ncol = 1)
  c2 <- matrix(stats::rnorm(2), nrow = 1, ncol = 2)
  constants <- c(2, as.numeric(c2 %*% c1))

  p <- Problem(Minimize(p1), list(a + p1 <= p2, b <= p3 + p3 + constants[1], c == constants[2]))

  # num_scalar variables
  n_variables <- size_metrics(p)@num_scalar_variables
  ref <- size(a) + size(b) + size(c)
  expect_equal(n_variables, ref)

  # num_scalar_data
  n_data <- size_metrics(p)@num_scalar_data
  # 2 and c2 %*% c1 are both single scalar constants
  ref <- prod(size(p1)) + prod(size(p2)) + prod(size(p3)) + length(constants)
  expect_equal(n_data, ref)

  # num_scalar_eq_constr
  n_eq_constr <- size_metrics(p)@num_scalar_eq_constr
  ref <- prod(dim(c2 %*% c1))
  expect_equal(n_eq_constr, ref)

  # num_scalar_leq_constr
  n_leq_constr <- size_metrics(p)@num_scalar_leq_constr
  ref <- prod(size(p3)) + prod(size(p2))
  expect_equal(n_leq_constr, ref)

  # max_data_dimension
  max_data_dim <- size_metrics(p)@max_data_dimension
  ref <- max(dim(p3))
  expect_equal(max_data_dim, ref)
})

test_that("Test the solver_stats method", {
  skip_on_cran()
  prob <- Problem(Minimize(p_norm(x)), list(x == 0))
  result <- solve(prob, solver = "ECOS")
  stats <- results$solver_stats
  expect_gt(stats@solve_time, 0)
  expect_gt(stats@setup_time, 0)
  expect_gt(stats@num_iters, 0)
  expect_true("info" %in% names(stats@extra_stats))
  
  prob <- Problem(Minimize(p_norm(x)), list(x == 0))
  result <- solve(prob, solver = "SCS")
  stats <- result$solver_stats
  expect_gt(stats@solve_time, 0)
  expect_gt(stats@setup_time, 0)
  expect_gt(stats@num_iters, 0)
  expect_true("info" %in% names(stats@extra_stats))
  
  prob <- Problem(Minimize(sum(x)), list(x == 0))
  result <- solve(prob, solver = "OSQP")
  stats <- result$solver_stats
  expect_gt(stats@solve_time, 0)
  
  # We do not populate setup_time for OSQP (OSQP decomposes time
  # into setup, solve, and polish; these are summed to obtain solve_time).
  expect_gt(stats@num_iters, 0)
  expect_true("info" %in% names(stats@extra_stats))
})

test_that("Test the get_problem_data method", {
  skip_on_cran()
  data <- get_problem_data(Problem(Minimize(exp(a) + 2)), "SCS")[[1]]
  dims <- data[[ConicSolver()@DIMS]]
  expect_equal(dims@exp, 1)
  expect_equal(length(data$c), 2)
  expect_equal(dim(data$A), c(3,2))

  data <- get_problem_data(Problem(Minimize(p_norm(x) + 3)), "ECOS")[[1]]
  dims <- data[[ConicSolver()@DIMS]]
  expect_equal(dims@soc, list(3))
  expect_equal(length(data$c), 3)
  expect_equal(dim(data$A), c(0,0))
  expect_equal(dim(data$G), c(3,3))

  if("CVXOPT" %in% INSTALLED_SOLVERS) {
    data <- get_problem_data(Problem(Minimize(p_norm(x) + 3)), "CVXOPT")[[1]]
    dims <- data[[ConicSolver()@DIMS]]
    expect_equal(dims@soc, list(3))
    # TODO: We cannot test whether the coefficients or offsets were correctly parsed until we update the CVXOPT interface.
  }
})

test_that("Test unpack results method", {
  skip_on_cran()
  prob <- Problem(Minimize(exp(a)), list(a == 0))
  tmp <- get_problem_data(prob, solver = "SCS")
  args <- tmp[[1]]
  chain <- tmp[[2]]
  inv <- tmp[[3]]
  data <- list(c = args$c, A = args$A, b = args$b)
  cones <- SCS.dims_to_solver_dict(args[[ConicSolver()@DIMS]])
  solution <- scs::scs(A = data$A, b = data$b, obj = data$c, cone = cones)
  prob <- Problem(Minimize(exp(a)), list(a == 0))
  result <- unpack_results(prob, solution, chain, inv)
  expect_equal(result$getValue(a), 0, tolerance = 1e-3)
  expect_equal(result$value, 1, tolerance = 1e-3)
  expect_equal(result$status, OPTIMAL)

  prob <- Problem(Minimize(p_norm(x)), list(x == 0))
  tmp <- get_problem_data(prob, solver = "ECOS")
  args <- tmp[[1]]
  chain <- tmp[[2]]
  inv <- tmp[[3]]
  cones <- ECOS.dims_to_solver_dict(args[[ConicSolver()@DIMS]])
  solution <- ECOSolveR::ECOS_csolve(args$c, args$G, args$h, cones, args$A, args$b)
  prob <- Problem(Minimize(p_norm(x)), list(x == 0))
  result <- unpack_results(prob, solution, chain, inv)
  expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = TOL)
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$status, OPTIMAL)
})

# test_that("Test silencing and enabling solver messages", {
# })

# test_that("Test registering other solve methods", {
#  Problem.register_solve("test", function(self) { 1 })
#  p <- Problem(Minimize(1))
#  result <- solve(p, method = "test")
#  expect_equal(result$value, 1)
#
#  test <- function(self, a, b = 2) { c(a, b) }
#  Problem.register_solve("test", test)
#  p <- Problem(Minimize(0))
#  result <- solve(p, 1, b = 3, method = "test")
#  expect_equal(result$value, c(1,3))
#  result <- solve(p, 1, method = "test")
#  expect_equal(result$value, c(1,2))
#  result <- solve(p, 1, method = "test", b = 4)
#  expect_equal(result$value, c(1,4))
# })

# TODO: Adapt this test to the reduction infrastructure.
# test_that("Test that variables and constraints keep a consistent order", {
#   num_solves <- 4
#   vars_lists <- list()
#   ineqs_lists <- list()
#   var_ids_order_created <- list()
#   for(k in 1:num_solves) {
#     sum <- 0
#     constraints <- list()
#     var_ids <- c()
#     for(i in 1:100) {
#       var <- Variable(name = as.character(i))
#       var_ids <- c(var_ids, var@id)
#       sum <- sum + var
#       constraints <- c(constraints, list(var >= i))
#     }
#     var_ids_order_created <- c(var_ids_order_created, list(var_ids))
#     obj <- Minimize(sum)
#     p <- Problem(obj, constraints)
#     canon <- canonicalize(p)
#     objective <- canon[[1]]
#     constraints <- canon[[2]]
#     sym_data <- SymData(objective, constraints, ECOS())
#
#     # Sort by offset
#     offsets <- sym_data@.var_offsets
#     vars_ <- as.numeric(names(sort(offsets, decreasing = FALSE)))
#     vars_lists <- c(vars_lists, list(vars_))
#     ineqs_lists <- c(ineqs_lists, list(sym_data@.constr_map[[LEQ_MAP]]))
#   }
#
#   # Verify order of variables is consistent
#   for(i in 1:num_solves)
#     expect_equal(var_ids_order_created[[i]], vars_lists[[i]])
#
#   for(i in 1:num_solves) {
#     idx <- 1
#     for(constr in ineqs_lists[[i]]) {
#       var_tmp <- get_expr_vars(constr$expr)[[1]]
#       var_id <- as.numeric(var_tmp[[1]])
#       expect_equal(var_ids_order_created[[i]][idx], var_id)
#       idx <- idx + 1
#     }
#   }
# })

# TODO: Adapt this test to the reduction infrastructure.
# test_that("Test removing duplicate constraints objects", {
#   eq <- (x == 2)
#   le <- (x <= 2)
#   obj <- 0
#
#   test <- function(self) {
#     canon <- canonicalize(self)
#     objective <- canon[[1]]
#     constraints <- canon[[2]]
#     sym_data <- SymData(objective, constraints)
#     list(length(sym_data@constr_map[EQ_MAP]), length(sym_data@constr_map[LEQ_MAP]))
#   }
#   # Problem.register_solve("test", test)
#   # p <- Problem(Minimize(obj), list(eq, eq, le, le))
#   # result <- solve(p, method = "test")
#   # expect_equal(result$value, c(1,1))
#
#   # Internal constraints
#   X <- Semidef(2)
#   obj <- sum_entries(X + X)
#   p <- Problem(Minimize(obj))
#   # result <- solve(p, method = "test")
#   # expect_equal(result$value, c(0,1))
#
#   # Duplicates from non-linear constraints
#   # exp <- norm(x, "2")
#   # prob <- Problem(Minimize(0), list(exp <= 1, exp <= 2))
#   # result <- solve(prob, method = "test")
#   # expect_equal(result$value, c(0,4))
# })

test_that("test the is_dcp method", {
  skip_on_cran()
  p <- Problem(Minimize(norm_inf(a)))
  expect_true(is_dcp(p))

  p <- Problem(Maximize(norm_inf(a)))
  expect_false(is_dcp(p))
  expect_error(solve(p))
})

test_that("test the is_qp method", {
  skip_on_cran()
  A <- matrix(rnorm(4*3), nrow = 4, ncol = 3)
  b <- matrix(rnorm(4), nrow = 4)
  Aeq <- matrix(rnorm(2*3), nrow = 2, ncol = 3)
  beq <- matrix(rnorm(2), nrow = 2, ncol = 1)
  Fmat <- matrix(rnorm(2*3), nrow = 2, ncol = 3)
  g <- matrix(rnorm(2), nrow = 2, ncol = 1)
  obj <- sum_squares(A %*% y - b)
  qpwa_obj <- 3*sum_squares(-abs(A %*% y)) + quad_over_lin(max_elemwise(abs(A %*% y), rep(3, 4)), 2)
  not_qpwa_obj <- 3*sum_squares(abs(A %*% y)) + quad_over_lin(min_elemwise(abs(A %*% y), rep(3, 4)), 2)

  p <- Problem(Minimize(obj), list())
  expect_true(is_qp(p))

  p <- Problem(Minimize(qpwa_obj), list())
  expect_true(is_qp(p))

  p <- Problem(Minimize(not_qpwa_obj), list())
  expect_false(is_qp(p))

  p <- Problem(Minimize(obj), list(Aeq %*% y == beq, Fmat %*% y <= g))
  expect_true(is_qp(p))

  p <- Problem(Minimize(obj), list(max_elemwise(1, 3*y) <= 200, abs(2*y) <= 100, p_norm(2*y, 1) <= 1000, Aeq %*% y == beq))
  expect_true(is_qp(p))

  p <- Problem(Minimize(qpwa_obj), list(max_elemwise(1, 3*y) <= 200, abs(2*y) <= 100, p_norm(2*y, 1) <= 1000, Aeq %*% y == beq))
  expect_true(is_qp(p))

  p <- Problem(Minimize(obj), list(max_elemwise(1, 3*y^2) <= 200))
  expect_false(is_qp(p))

  p <- Problem(Minimize(qpwa_obj), list(max_elemwise(1, 3*y^2) <= 200))
  expect_false(is_qp(p))
})

test_that("test problems involving variables with the same name", {
  skip_on_cran()
  var <- Variable(name = "a")
  p <- Problem(Maximize(a + var), list(var == 2 + a, var <= 3))
  result <- solve(p)
  expect_equal(result$value, 4.0, tolerance = TOL)
  expect_equal(result$getValue(a), 1, tolerance = TOL)
  expect_equal(result$getValue(var), 3, tolerance = TOL)
})

test_that("test adding problems", {
  skip_on_cran()
  prob1 <- Problem(Minimize(a), list(a >= b))
  prob2 <- Problem(Minimize(2*b), list(a >= 1, b >= 2))
  prob_minimize <- prob1 + prob2
  expect_equal(length(prob_minimize@constraints), 3)
  result <- solve(prob_minimize)
  expect_equal(result$value, 6, tolerance = TOL)

  prob3 <- Problem(Maximize(a), list(b <= 1))
  prob4 <- Problem(Maximize(2*b), list(a <= 2))
  prob_maximize <- prob3 + prob4
  expect_equal(length(prob_maximize@constraints), 2)
  result <- solve(prob_maximize)
  expect_equal(result$value, 4, tolerance = TOL)

  # Test using the sum function
  prob5 <- Problem(Minimize(3*a))
  prob_sum <- Reduce("+", list(prob1, prob2, prob5))
  expect_equal(length(prob_sum@constraints), 3)
  result <- solve(prob_sum)
  expect_equal(result$value, 12, tolerance = TOL)
  prob_sum <- Reduce("+", list(prob1))
  expect_equal(length(prob_sum@constraints), 1)

  # Test Minimize + Maximize
  expect_error(prob_bad_sum <- prob1 + prob3, "Problem does not follow DCP rules")
})

test_that("test problem multiplication by scalar", {
  skip_on_cran()
  prob1 <- Problem(Minimize(a^2), list(a >= 2))
  answer <- solve(prob1)$value
  factors <- c(0, 1, 2.3, -4.321)
  for(f in factors) {
    expect_equal(solve(f * prob1)$value, f * answer, tolerance = TOL)
    expect_equal(solve(prob1 * f)$value, f * answer, tolerance = TOL)
  }
})

test_that("test problem linear combinations", {
  skip_on_cran()
  prob1 <- Problem(Minimize(a), list(a >= b))
  prob2 <- Problem(Minimize(2*b), list(a >= 1, b >= 2))
  prob3 <- Problem(Maximize(-(b + a)^2), list(b >= 3))

  # Simple addition and multiplication.
  combo1 <- prob1 + 2 * prob2
  combo1_ref <- Problem(Minimize(a + 4*b), list(a >= b, a >= 1, b >= 2))
  expect_equal(solve(combo1)$value, solve(combo1_ref)$value)

  # Division and subtraction.
  combo2 <- prob1 - prob3/2
  combo2_ref <- Problem(Minimize(a + (b + a)^2/2), list(b >= 3, a >= b))
  expect_equal(solve(combo2)$value, solve(combo2_ref)$value)

  # Multiplication with 0 (prob2's constraints should still hold).
  combo3 <- prob1 + 0*prob2 - 3*prob3
  combo3_ref <- Problem(Minimize(a + 3*(b + a)^2), list(a >= b, a >= 1, b >= 3))
  expect_equal(solve(combo3)$value, solve(combo3_ref)$value, tolerance = TOL)
})

# test_that("test solving problems in parallel", {
#   p <- Parameter()
#   problem <- Problem(Minimize(a^2 + b^2 + p), list(b >= 2, a >= 1))
#   value(p) <- 1
#
#   # Ensure that parallel solver still works after repeated calls
#   for(i in 1:2) {
#     result <- solve(problem, parallel = TRUE)
#     expect_equal(result$value, 6.0, tolerance = TOL)
#     expect_equal(result$status, "optimal")
#     expect_equal(result$getValue(a), 1, tolerance = TOL)
#     expect_equal(result$getValue(b), 2, tolerance = TOL)
#   }
#
#   # The constant p should not be a separate problem, but rather added to the first separable problem
#   expect_true(length(problem@.separable_problems) == 2)
#
#   # Ensure that parallel solver works with options
#   result <- solve(problem, parallel = TRUE, verbose = TRUE, warm_start = TRUE)
#   expect_equal(result$value, 6.0, tolerance = TOL)
#   expect_equal(result$status, "optimal")
#   expect_equal(result$getValue(a), 1, tolerance = TOL)
#   expect_equal(result$getValue(b), 2, tolerance = TOL)
#
#   # Ensure that parallel solver works when problem changes
#   objective <- Minimize(a^2 + b^2)
#   problem <- Problem(objective, problem@constraints)
#   result <- solve(problem, parallel = TRUE)
#   expect_equal(result$value, 5.0, tolerance = TOL)
#   expect_equal(result$status, "optimal")
#   expect_equal(result$getValue(a), 1, tolerance = TOL)
#   expect_equal(result$getValue(b), 2, tolerance = TOL)
# })

test_that("Test scalar LP problems", {
  skip_on_cran()
  p <- Problem(Minimize(3*a), list(a >= 2))
  result <- solve(p)
  expect_equal(result$value, 6, tolerance = TOL)
  expect_equal(result$getValue(a), 2, tolerance = TOL)

  p <- Problem(Maximize(3*a - b), list(a <= 2, b == a, b <= 5))
  result <- solve(p)
  expect_equal(result$value, 4.0, tolerance = TOL)
  expect_equal(result$getValue(a), 2, tolerance = TOL)
  expect_equal(result$getValue(b), 2, tolerance = TOL)

  # With a constant in the objective
  p <- Problem(Minimize(3*a - b + 100), list(a >= 2, b + 5*c - 2 == a, b <= 5 + c))
  result <- solve(p)
  expect_equal(result$value, 101+1.0/6, tolerance = TOL)
  expect_equal(result$getValue(a), 2, tolerance = TOL)
  expect_equal(result$getValue(b), 5-1.0/6, tolerance = TOL)
  expect_equal(result$getValue(c), -1.0/6, tolerance = TOL)

  # Test status and value
  exp <- Maximize(a)
  p <- Problem(exp, list(a <= 2))
  result <- solve(p, solver = "ECOS")
  # expect_equal(result$value, value(p))
  expect_equal(result$status, "optimal")
  expect_false(is.na(result$getValue(a)))
  expect_false(is.na(result$getDualValue(p@constraints[[1]])))

  # Unbounded problems
  p <- Problem(Maximize(a), list(a >= 2))
  result <- solve(p, solver = "ECOS")
  expect_equal(result$value, Inf)
  expect_true(is.na(result$getValue(a)))
  expect_true(is.na(result$getDualValue(p@constraints[[1]])))

  if("CVXOPT" %in% INSTALLED_SOLVERS) {
    p <- Problem(Minimize(-a), list(a >= 2))
    result <- solve(p, solver = "CVXOPT")
    # expect_equal(result$value, value(p))
    expect_equal(result$status, "unbounded")
    expect_equal(result$value, -Inf)
  }

  # Infeasible problems
  p <- Problem(Maximize(a), list(a >= 2, a <= 1))
  # a <- save_value(a, 2)
  # p@constraints[[1]] <- save_value(p@constraints[[1]], 2)

  result <- solve(p, solver = "ECOS")
  # expect_equal(result$value, value(p))
  expect_equal(result$status, "infeasible")
  expect_equal(result$value, -Inf)
  expect_true(is.na(result$getValue(a)))
  expect_true(is.na(result$getDualValue(p@constraints[[1]])))

  p <- Problem(Minimize(-a), list(a >= 2, a <= 1))
  result <- solve(p, solver = "ECOS")
  # expect_equal(result$value, value(p))
  expect_equal(result$status, "infeasible")
  expect_equal(result$value, Inf)
})

test_that("Test vector LP problems", {
  skip_on_cran()
  c <- Constant(matrix(c(1, 2), nrow = 2, ncol = 1))@value
  p <- Problem(Minimize(t(c) %*% x), list(x >= c))
  result <- solve(p)
  expect_equal(result$value, 5, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(1, 2)), tolerance = TOL)

  A <- Constant(rbind(c(3, 5), c(1, 2)))@value
  I <- Constant(rbind(c(1, 0), c(0, 1)))
  p <- Problem(Minimize(t(c) %*% x + a), list(A %*% x >= c(-1, 1), 4*I %*% z == x, z >= c(2,2), a >= 2))
  result <- solve(p)
  expect_equal(result$value, 26, tolerance = 1e-3)
  obj <- result$getValue(t(c) %*% x + a)
  expect_equal(obj, result$value, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(8,8)), tolerance = 1e-3)
  expect_equal(result$getValue(z), matrix(c(2,2)), tolerance = 1e-3)
})

test_that("Test ECOS with no inequality constraints", {
  skip_on_cran()
  Tmat <- value(Constant(matrix(1, nrow = 2, ncol = 2)))
  p <- Problem(Minimize(1), list(A == Tmat))
  result <- solve(p, solver = "ECOS")
  expect_equal(result$value, 1, tolerance = TOL)
  expect_equal(result$getValue(A), Tmat, tolerance = TOL)
})

test_that("Test matrix LP problems", {
  skip_on_cran()
  Tmat <- value(Constant(matrix(1, nrow = 2, ncol = 2)))
  p <- Problem(Minimize(1), list(A == Tmat))
  result <- solve(p)
  expect_equal(result$value, 1, tolerance = TOL)
  expect_equal(result$getValue(A), Tmat, tolerance = TOL)

  Tmat <- value(Constant(matrix(1, nrow = 2, ncol = 3)*2))
  p <- Problem(Minimize(1), list(A >= Tmat %*% C, A == B, C == t(Tmat)))
  result <- solve(p)
  expect_equal(result$value, 1, tolerance = TOL)
  expect_equal(result$getValue(A), result$getValue(B), tolerance = TOL)
  expect_equal(result$getValue(C), t(Tmat), tolerance = TOL)
  expect_true(all(result$getValue(A) >= Tmat %*% result$getValue(C)))

  # Test variables are dense
  expect_true(is(result$getValue(A), "matrix"))
})

test_that("Test variable promotion", {
  skip_on_cran()
  p <- Problem(Minimize(a), list(x <= a, x == c(1, 2)))
  result <- solve(p)
  expect_equal(result$value, 2, tolerance = TOL)
  expect_equal(result$getValue(a), 2, tolerance = TOL)

  p <- Problem(Minimize(a), list(A <= a, A == rbind(c(1,2), c(3,4))))
  result <- solve(p)
  expect_equal(result$value, 4, tolerance = TOL)
  expect_equal(result$getValue(a), 4, tolerance = TOL)

  # Promotion must happen before multiplication
  p <- Problem(Minimize(matrix(1, nrow = 1, ncol = 2) %*% (x + a + 1)), list(a + x >= c(1, 2)))
  result <- solve(p)
  expect_equal(result$value, 5, tolerance = TOL)
})

test_that("Test parameter promotion", {
  skip_on_cran()
  a <- Parameter()
  value(a) <- 2
  exp <- cbind(c(1,2), c(3,4))*a
  expect_false(any(value(exp) - 2*cbind(c(1,2), c(3,4)) != 0))
})

test_that("test problems with parameters", {
  skip_on_cran()
  p1 <- Parameter()
  p2 <- Parameter(3, nonpos = TRUE)
  p3 <- Parameter(4, 4, nonneg = TRUE)
  p <- Problem(Maximize(p1*a), list(a + p1 <= p2, b <= p3 + p3 + 2))

  value(p1) <- 2
  value(p2) <- -matrix(1, nrow = 3, ncol = 1)
  value(p3) <- matrix(1, nrow = 4, ncol = 4)
  ##p <- Problem(Maximize(p1*a), list(a + p1 <= p2, b <= p3 + p3 + 2))
  result <- solve(p)
  expect_equal(result$value, -6, tolerance = TOL)

  value(p1) <- NA_real_
  p <- Problem(Maximize(p1*a), list(a + p1 <= p2, b <= p3 + p3 + 2))
  expect_error(solve(p))
})

test_that("test problems with norm_inf", {
  skip_on_cran()
  # Constant argument
  p <- Problem(Minimize(norm_inf(-2)))
  result <- solve(p)
  expect_equal(result$value, 2, tolerance = TOL)

  # Scalar arguments
  p <- Problem(Minimize(norm_inf(a)), list(a >= 2))
  result <- solve(p)
  expect_equal(result$value, 2, tolerance = TOL)
  expect_equal(result$getValue(a), 2, tolerance = TOL)

  p <- Problem(Minimize(3*norm_inf(a + 2*b) + c), list(a >= 2, b <= -1, c == 3))
  result <- solve(p)
  expect_equal(result$value, 3, tolerance = TOL)
  expect_equal(result$getValue(a + 2*b), 0, tolerance = TOL)
  expect_equal(result$getValue(c), 3, tolerance = TOL)

  # Maximize
  p <- Problem(Maximize(-norm_inf(a)), list(a <= -2))
  result <- solve(p)
  expect_equal(result$value, -2, tolerance = TOL)
  expect_equal(result$getValue(a), -2, tolerance = TOL)

  # Vector arguments
  p <- Problem(Minimize(norm_inf(x - z) + 5), list(x >= c(2,3), z <= c(-1,-4)))
  result <- solve(p, solver = "ECOS")
  expect_equal(result$value, 12, tolerance = TOL)
  expect_equal(result$getValue(x[2] - z[2]), 7, tolerance = TOL)
})

test_that("Test problems with norm1", {
  skip_on_cran()
  # Constant argument
  p <- Problem(Minimize(norm1(-2)))
  result <- solve(p)
  expect_equal(result$value, 2, tolerance = TOL)

  # Scalar arguments
  p <- Problem(Minimize(norm1(a)), list(a <= -2))
  result <- solve(p)
  expect_equal(result$value, 2, tolerance = TOL)
  expect_equal(result$getValue(a), -2, tolerance = TOL)

  # Maximize
  p <- Problem(Maximize(-norm1(a)), list(a <= -2))
  result <- solve(p)
  expect_equal(result$value, -2, tolerance = TOL)
  expect_equal(result$getValue(a), -2, tolerance = TOL)

  # Vector arguments
  p <- Problem(Minimize(norm1(x - z) + 5), list(x >= c(2,3), z <= c(-1,-4)))
  result <- solve(p)
  expect_equal(result$value, 15, tolerance = TOL)
  expect_equal(result$getValue(x[2] - z[2]), 7, tolerance = TOL)
})

test_that("Test problems with norm2", {
  skip_on_cran()
  # Constant argument
  p <- Problem(Minimize(norm2(-2)))
  result <- solve(p)
  expect_equal(result$value, 2, tolerance = TOL)

  # Scalar arguments
  p <- Problem(Minimize(norm2(a)), list(a <= -2))
  result <- solve(p)
  expect_equal(result$value, 2, tolerance = TOL)
  expect_equal(result$getValue(a), -2, tolerance = TOL)

  # Maximize
  p <- Problem(Maximize(-norm2(a)), list(a <= -2))
  result <- solve(p)
  expect_equal(result$value, -2, tolerance = TOL)
  expect_equal(result$getValue(a), -2, tolerance = TOL)

  # Vector arguments
  p <- Problem(Minimize(norm2(x - z) + 5), list(x >= c(2,3), z <= c(-1,-4)))
  result <- solve(p)
  expect_equal(result$value, 12.61577, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(2,3)), tolerance = TOL)
  expect_equal(result$getValue(z), matrix(c(-1,-4)), tolerance = TOL)

  # Row arguments
  p <- Problem(Minimize(norm2(t(x - z)) + 5), list(x >= c(2,3), z <= c(-1,-4)))
  result <- solve(p)
  expect_equal(result$value, 12.61577, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(2,3)), tolerance = TOL)
  expect_equal(result$getValue(z), matrix(c(-1,-4)), tolerance = TOL)
})

test_that("Test problems with abs", {
  skip_on_cran()
  p <- Problem(Minimize(sum_entries(abs(A))), list(-2 >= A))
  result <- solve(p)
  expect_equal(result$value, 8, tolerance = TOL)
  expect_equal(result$getValue(A), matrix(rep(-2, 4), nrow = 2, ncol = 2), tolerance = TOL)
})

test_that("Test problems with quad_form", {
  skip_on_cran()
  expect_error(solve(Problem(Minimize(quad_form(x, A)))), "At least one argument to QuadForm must be constant.")
  expect_error(solve(Problem(Minimize(quad_form(1, A)))), "Invalid dimensions for arguments.")
  expect_error(solve(Problem(Minimize(quad_form(x, rbind(c(-1,0), c(0,9)))))))              # Error: Problem does not follow DCP rules

  P <- rbind(c(4,0), c(0,9))
  p <- Problem(Minimize(quad_form(x, P)), list(x >= 1))
  result <- solve(p)
  expect_equal(result$value, 13, tolerance = 1e-3)

  c <- c(1,2)
  p <- Problem(Minimize(quad_form(c, A)), list(A >= 1))
  result <- solve(p)
  expect_equal(result$value, 9, tolerance = TOL)

  c <- c(1,2)
  P <- rbind(c(4,0), c(0,9))
  p <- Problem(Minimize(quad_form(c, P)))
  result <- solve(p)
  expect_equal(result$value, 40)
})

test_that("Test combining atoms", {
  skip_on_cran()
  p <- Problem(Minimize(norm2(5 + p_norm(z,1) + p_norm(x,1) + norm_inf(x - z))), list(x >= c(2,3), z <= c(-1,-4), p_norm(x + z,2) <= 2))
  result <- solve(p)
  expect_equal(result$value, 22, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(2,3)), tolerance = TOL)
  expect_equal(result$getValue(z), matrix(c(-1,-4)), tolerance = TOL)
})

test_that("Test multiplying by constant atoms", {
  skip_on_cran()
  p <- Problem(Minimize(norm2(c(3,4)) * a), list(a >= 2))
  result <- solve(p)
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equal(result$getValue(a), 2, tolerance = TOL)
})

test_that("Test recovery of dual variables", {
  skip_on_cran()
  for(solver in c("ECOS", "SCS", "CVXOPT")) {
    if(solver %in% INSTALLED_SOLVERS) {
      if(solver == "SCS")
        acc <- 1e-1
      else
        acc <- 1e-5
      p <- Problem(Minimize(p_norm(x + z, 1)), list(x >= c(2,3), cbind(c(1,2), c(3,4)) %*% z == c(-1,-4), p_norm(x + z, 2) <= 100))
      result <- solve(p, solver = solver)
      expect_equal(result$value, 4, tolerance = acc)
      expect_equal(result$getValue(x), matrix(c(4,3)), tolerance = acc)
      expect_equal(result$getValue(z), matrix(c(-4,1)), tolerance = acc)

      # Dual values
      expect_equal(result$getDualValue(p@constraints[[1]]), matrix(c(0,1)), tolerance = acc)
      expect_equal(result$getDualValue(p@constraints[[2]]), matrix(c(-1,0.5)), tolerance = acc)
      expect_equal(result$getDualValue(p@constraints[[3]]), 0, tolerance = acc)

      Tmat <- matrix(1, nrow = 2, ncol = 3) * 2
      c <- matrix(c(3,4), nrow = 1, ncol = 2)
      p <- Problem(Minimize(1), list(A >= Tmat %*% C, A == B, C == t(Tmat)))
      result <- solve(p, solver = solver)

      # Dual values
      expect_equal(result$getDualValue(p@constraints[[1]]), matrix(0, nrow = 2, ncol = 2), tolerance = acc)
      expect_equal(result$getDualValue(p@constraints[[2]]), matrix(0, nrow = 2, ncol = 2), tolerance = acc)
      expect_equal(result$getDualValue(p@constraints[[3]]), matrix(0, nrow = 3, ncol = 2), tolerance = acc)
    }
  }
})

test_that("Test problems with indexing", {
  skip_on_cran()
  
  # Vector variables.
  p <- Problem(Maximize(x[1]), list(x[1] <= 2, x[2] == 3))
  result <- solve(p, solver = "SCS", eps = 1e-8)
  expect_equal(result$value, 2, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(2,3)), tolerance = TOL)

  n <- 10
  Aloc <- matrix(0:(n^2-1), nrow = n, ncol = n, byrow = TRUE)
  xloc <- Variable(n,n)
  p <- Problem(Minimize(sum_entries(xloc)), list(xloc == Aloc))
  result <- solve(p, solver = "SCS", eps = 1e-8)
  answer <- n*n*(n*n+1)/2 - n*n
  expect_equal(result$value, answer, tolerance = TOL)

  # Matrix variables.
  obj <- A[1,1] + A[1,2] + A[2,2] + A[2,1]
  p <- Problem(Maximize(obj), list(A <= rbind(c(1,-2), c(-3,4))))
  result <- solve(p, solver = "SCS", eps = 1e-8)
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(A), rbind(c(1,-2), c(-3,4)))

  # Indexing arithmetic expressions.
  exp <- cbind(c(1,2), c(3,4)) %*% z + x
  p <- Problem(Minimize(exp[2]), list(x == z, z == c(1,2)))
  result <- solve(p, solver = "SCS", eps = 1e-8)
  expect_equal(result$value, 12, tolerance = TOL)
  expect_equal(result$getValue(x), result$getValue(z), tolerance = TOL)
})

test_that("Test problems that have special types as indices", {
  # Test with default 32-bit int indices.
  cost <- x[1:integer(2)][1]
  p <- Problem(Minimize(cost), list(x == 1))
  result <- solve(p, solver = "SCS")
  expect_equal(result$value, 1, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(1,1)), tolerance = TOL)
  
  # Test with 64-bit int indices.
  if(require(bit64)) {
    cost <- x[1:integer64(2)][1]
    p <- Problem(Minimize(cost), list(x == 1))
    result <- solve(p, solver = "SCS")
    expect_equal(result$value, 1, tolerance = TOL)
    expect_equal(result$getValue(x), matrix(c(1,1)), tolerance = TOL)
  } else {
    print("bit64 library not installed. Skipping test.")
    return()
  }
})

test_that("Test problems with slicing", {
  skip_on_cran()
  p <- Problem(Maximize(sum_entries(C)), list(C[2:3,] <= 2, C[1,] == 1))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equal(result$getValue(C), rbind(c(1,2,2), c(1,2,2)))

  p <- Problem(Maximize(sum_entries(C[seq(1,3,2),2])), list(C[2:3,] <= 2, C[1,] == 1))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 3, tolerance = TOL)
  expect_equal(result$getValue(C[seq(1,3,2),2]), matrix(c(1,2)), tolerance = TOL)

  p <- Problem(Maximize(sum_entries((C[1:2,] + A)[,1:2])), 
               list(C[2:3,] <= 2, C[1,] == 1, (A + B)[,1] == 3,
                    (A + B)[,2] == 2, B == 1))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 12, tolerance = TOL)
  expect_equal(result$getValue(C[1:2,]), rbind(c(1,2), c(1,2)), tolerance = TOL)
  expect_equal(result$getValue(A), rbind(c(2,2), c(1,1)), tolerance = TOL)

  p <- Problem(Maximize(matrix(c(3,4), nrow = 1, ncol = 2) %*% (C[1:2,] + A)[,1]),
               list(C[2:3,] <= 2, C[1,] == 1, matrix(c(1,2), nrow = 1, ncol = 2) %*% (A + B)[,1] == 3,
                    (A + B)[,2] == 2, B == 1, 3*A[,1] <= 3))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 12, tolerance = TOL)
  expect_equal(result$getValue(C[1:2,1]), matrix(c(1,2)), tolerance = TOL)
  expect_equal(result$getValue(A), rbind(c(1,-0.5), c(1,1)), tolerance = TOL)

  p <- Problem(Minimize(p_norm((C[1:2,] + A)[,1]), p = 2), 
               list(C[2:3,] <= 2, C[1,] == 1, (A + B)[,1] == 3, (A + B)[,2] == 2, B == 1))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 3, tolerance = TOL)
  expect_equal(result$getValue(C[1:2,1]), matrix(c(1,-2)), tolerance = 1e-3)
  expect_equal(result$getValue(A), rbind(c(2,2), c(1,1)), tolerance = TOL)

  # Transpose of slice.
  p <- Problem(Maximize(sum_entries(C)), list(t(C[2:3,]) <= 2, t(C[1,]) == 1))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equal(result$getValue(C), rbind(c(1,2,2), c(1,2,2)), tolerance = TOL)
})

test_that("Test the vstack function", {
  skip_on_cran()
  a <- Variable(1, 1, name = "a")
  b <- Variable(1, 1, name = "b")
  
  x <- Variable(2, 1, name = "x")
  y <- Variable(3, 1, name = "y")
  
  c <- matrix(1, nrow = 1, ncol = 5)
  p <- Problem(Minimize(c %*% vstack(x, y)), list(x == c(1,2), y == c(3,4,5)))
  result <- solve(p, solver = "SCS", eps = 1e-5)
  expect_equal(result$value, 15, tolerance = TOL)

  c <- matrix(1, nrow = 1, ncol = 4)
  p <- Problem(Minimize(c %*% vstack(x, x)), list(x == c(1,2)))
  result <- solve(p, solver = "SCS", eps = 1e-8)
  expect_equal(result$value, 6, tolerance = TOL)

  c <- matrix(1, nrow = 2, ncol = 2)
  p <- Problem(Minimize(sum_entries(vstack(A, C))), list(A >= 2*c, C == -2))
  result <- solve(p, solver = "SCS", eps = 1e-8)
  expect_equal(result$value, -4, tolerance = TOL)

  c <- matrix(1, nrow = 1, ncol = 2)
  p <- Problem(Minimize(sum_entries(vstack(c %*% A, c %*% B))), list(A >= 2, B == -2))
  result <- solve(p, solver = "SCS", eps = 1e-8)
  expect_equal(result$value, 0, tolerance = TOL)

  c <- matrix(c(1,-1), nrow = 2, ncol = 1)
  p <- Problem(Minimize(t(c) %*% vstack(square(a), sqrt(b))), list(a == 2, b == 16))
  expect_error(solve(p, solver = "SCS", eps = 1e-5), "Problem does not follow DCP rules.")
})

test_that("Test the hstack atom", {
  skip_on_cran()
  a <- Variable(1, 1, name = "a")
  b <- Variable(1, 1, name = "b")
  
  x <- Variable(2, 1, name = "x")
  y <- Variable(3, 1, name = "y")
  
  c <- matrix(1, nrow = 1, ncol = 5)
  p <- Problem(Minimize(c %*% t(hstack(t(x), t(y)))), list(x == c(1,2), y == c(3,4,5)))
  result <- solve(p, solver = "SCS", eps = 1e-8)
  expect_equal(result$value, 15, tolerance = TOL)

  c <- matrix(1, nrow = 1, ncol = 4)
  p <- Problem(Minimize(c %*% t(hstack(t(x), t(x)))), list(x == c(1,2)))
  result <- solve(p, solver = "SCS", eps = 1e-8)
  expect_equal(result$value, 6, tolerance = TOL)

  c <- matrix(1, nrow = 2, ncol = 2)
  p <- Problem(Minimize(sum_entries(hstack(t(A), t(C)))), list(A >= 2*c, C == -2))
  result <- solve(p, solver = "SCS", eps = 1e-8)
  expect_equal(result$value, -4, tolerance = TOL)

  D <- Variable(3,3)
  expr <- hstack(C, D)
  p <- Problem(Minimize(expr[1,2] + sum_entries(hstack(expr, expr))), list(C >= 0, D >= 0, D[1,1] == 2, C[1,2] == 3))
  result <- solve(p, solver = "SCS", eps = 1e-8)
  expect_equal(result$value, 13, tolerance = TOL)

  c <- matrix(c(1,-1), nrow = 2, ncol = 1)
  p <- Problem(Minimize(t(c) %*% t(hstack(t(square(a)), t(sqrt(b))))), list(a == 2, b == 16))
  expect_error(solve(p, solver = "SCS", eps = 1e-5), "Problem does not follow DCP rules.")
})

test_that("Test using a CVXR expression as an objective", {
  skip_on_cran()
  expect_error(Problem(x+2), "Problem objective must be Minimize or Maximize.", fixed = TRUE)
})

test_that("Test variable transpose", {
  skip_on_cran()
  p <- Problem(Minimize(sum_entries(x)), list(t(x) >= matrix(c(1,2), nrow = 1, ncol = 2)))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 3, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(1,2)), tolerance = TOL)

  p <- Problem(Minimize(sum_entries(C)), list(matrix(c(1,1), nrow = 1, ncol = 2) %*% t(C) >= matrix(0:2, nrow = 1, ncol = 3)))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  value <- result$getValue(C)

  constraints <- lapply(1:3, function(i) { 1*C[i,1] + 1*C[i,2] >= (i-1) })
  p <- Problem(Minimize(sum_entries(C)), constraints)
  result2 <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, result2$value, tolerance = TOL)
  expect_equal(result2$getValue(C), value, tolerance = TOL)

  p <- Problem(Minimize(A[1,2] - t(A)[2,1]), list(A == rbind(c(1,2), c(3,4))))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 0, tolerance = TOL)
  
  p <- Problem(Minimize(sum_entries(x)), list(t(-x) <= 1))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, -2, tolerance = TOL)

  c <- matrix(c(1,-1), nrow = 2, ncol = 1)
  p <- Problem(Minimize(max_elemwise(t(c), 2, 2 + t(c))[1,2]))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 2, tolerance = TOL)

  c <- cbind(c(1,-1,2), c(1,-1,2))
  p <- Problem(Minimize(sum_entries(t(max_elemwise(c, 2, 2+c))[,1])))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 6, tolerance = TOL)

  c <- cbind(c(1,-1,2), c(1,-1,2))
  p <- Problem(Minimize(sum_entries(t(square(t(c)))[,1])))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 6, tolerance = TOL)

  # Slice of transpose.
  p <- Problem(Maximize(sum_entries(C)), list(t(C)[,2:3] <= 2, t(C)[,1] == 1))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equal(result$getValue(C), cbind(c(1,2,2), c(1,2,2)), tolerance = TOL)
})

test_that("Test multiplication on the left by a non-constant", {
  skip_on_cran()
  c <- matrix(c(1,2), nrow = 2, ncol = 1)
  p <- Problem(Minimize(t(c) %*% A %*% c), list(A >= 2))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 18, tolerance = TOL)

  p <- Problem(Minimize(a*2), list(a >= 2))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 4, tolerance = TOL)

  p <- Problem(Minimize(t(x) %*% c), list(x >= 2))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 6, tolerance = TOL)

  p <- Problem(Minimize((t(x) + t(z)) %*% c), list(x >= 2, z >= 1))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 9, tolerance = TOL)
  
  A <- matrix(1, nrow = 5, ncol = 10)
  x <- Variable(5)
  p <- Problem(Minimize(sum_entries(t(x) %*% A)), list(x >= 0))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 0, tolerance = TOL)
})

test_that("Test redundant constraints", {
  skip_on_cran()
  obj <- Minimize(sum_entries(x))
  constraints <- list(x == 2, x == 2, t(x) == 2, x[1] == 2)
  p <- Problem(obj, constraints)
  result <- solve(p, solver = "SCS")
  expect_equal(result$value, 4, tolerance = TOL)

  obj <- Minimize(sum_entries(square(x)))
  constraints <- list(x == x)
  p <- Problem(obj, constraints)
  result <- solve(p, solver = "SCS")
  expect_equal(result$value, 0, tolerance = TOL)
  
  obj <- Minimize(sum(square(x)))
  constraints <- list(x == x)
  problem <- Problem(obj, constraints)
  expect_error(solve(problem, solver = "ECOS"), 
               "ECOS cannot handle sparse data with nnz == 0; this is a bug in ECOS, and it indicates that your problem might have redundant constraints.", fixed = TRUE)
})

test_that("Test that SDP symmetry is enforced", {
  skip_on_cran()
  p <- Problem(Minimize(lambda_max(A)), list(A >= 2))
  result <- solve(p, solver = "SCS")
  expect_equal(result$getValue(A), t(result$getValue(A)), tolerance = 1e-3)

  p <- Problem(Minimize(lambda_max(A)), list(A == rbind(c(1,2), c(3,4))))
  result <- solve(p, solver = "SCS")
  expect_equal(result$status, INFEASIBLE)
})

test_that("Test SDP", {
  skip_on_cran()
  
  # Ensure SDP constraints enforce transpose.
  obj <- Maximize(A[2,1] - A[1,2])
  p <- Problem(obj, list(lambda_max(A) <= 100, A[1,1] == 2, A[2,2] == 2, A[2,1] == 2))
  result <- solve(p, solver = "SCS")
  expect_equal(result$value, 0, tolerance = 1e-3)
})

test_that("Test getting values for expressions", {
  skip_on_cran()
  diff_exp <- x - z
  inf_exp <- norm_inf(diff_exp)
  sum_entries_exp <- 5 + norm1(z) + norm1(x) + inf_exp
  constr_exp <- p_norm(x + z, p = 2)
  obj <- p_norm(sum_entries_exp, p = 2)
  p <- Problem(Minimize(obj), list(x >= c(2,3), z <= c(-1,-4), constr_exp <= 2))
  result <- solve(p, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 22, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(2,3)), tolerance = TOL)
  expect_equal(result$getValue(z), matrix(c(-1,-4)), tolerance = TOL)

  # Expression values.
  xs <- result$getValue(x)
  zs <- result$getValue(z)
  expect_equal(result$getValue(diff_exp), xs - zs, tolerance = TOL)
  expect_equal(result$getValue(inf_exp), norm(xs - zs, "I"))
  expect_equal(result$getValue(sum_entries_exp), 5 + norm(zs, "1") + norm(xs, "1") + norm(xs - zs, "I"))
  expect_equal(result$getValue(constr_exp), norm(xs + zs, "2"))
  expect_equal(result$getValue(obj), result$value)
})

test_that("Test multiplication by zero", {
  skip_on_cran()
  value(a) <- 1
  expr <- 0*a
  expect_equal(value(expr), matrix(0))
  obj <- Minimize(expr)
  p <- Problem(obj)
  result <- solve(p, solver = "ECOS")
  expect_equal(result$value, 0, tolerance = TOL)
  expect_false(is.na(result$getValue(a)))
})

test_that("Tests a problem with division", {
  skip_on_cran()
  obj <- Minimize(norm_inf(A/5))
  p <- Problem(obj, list(A >= 5))
  result <- solve(p, solver = "ECOS")
  expect_equal(result$value, 1, tolerance = TOL)
  
  c <- Constant(rbind(c(1,-1), c(2,-2)))
  expr <- A/(1/c)
  obj <- Minimize(norm_inf(expr))
  p <- Problem(obj, list(A == 5))
  result <- solve(p, solver = "ECOS")
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equal(result$getValue(expr), rbind(c(5,-5), c(10,-10)), tolerance = TOL)
  
  # Test with a sparse matrix.
  c <- Matrix(c(1,2), sparse = TRUE)
  c <- Constant(c)
  expr <- Minimize(norm_inf(expr))
  p <- Problem(obj, list(x == 5))
  result <- solve(p, solver = "ECOS")
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equivalent(as.matrix(result$getValue(expr)), matrix(c(5,10)), tolerance = TOL)
  
  # Test promotion.
  c <- rbind(c(1,-1), c(2,-2))
  c <- Constant(c)
  expr <- a/(1/c)
  obj <- Minimize(norm_inf(expr))
  p <- Problem(obj, list(a == 5))
  result <- solve(p, solver = "ECOS")
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equivalent(result$getValue(expr), rbind(c(5,-5), c(10,-10)), tolerance = TOL)
})

test_that("Tests problems with multiply", {
  skip_on_cran()
  c <- rbind(c(1,-1), c(2,-2))
  expr <- multiply(c, A)
  obj <- Minimize(norm_inf(expr))
  p <- Problem(obj, list(A == 5))
  result <- solve(p, solver = "ECOS")
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equal(result$getValue(expr), rbind(c(5,-5), c(10,-10)), tolerance = TOL)

  # Test with a sparse matrix.
  c <- Matrix(c(1,2), sparse = TRUE)
  expr <- multiply(c, x)
  obj <- Minimize(norm_inf(expr))
  p <- Problem(obj, list(x == 5))
  result <- solve(p, solver = "ECOS")
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equivalent(as.matrix(result$getValue(expr)), matrix(c(5,10)), tolerance = TOL)

  # Test promotion.
  c <- rbind(c(1,-1), c(2,-2))
  expr <- multiply(c, a)
  obj <- Minimize(norm_inf(expr))
  p <- Problem(obj, list(a == 5))
  result <- solve(p, solver = "SCS")
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equivalent(as.matrix(result$getValue(expr)), rbind(c(5,-5), c(10,-10)), tolerance = TOL)
})

test_that("Tests that errors occur when you use an invalid solver", {
  skip_on_cran()
  expect_error(solve(Problem(Minimize(Variable(boolean = TRUE))), solver = "ECOS"))
  expect_error(solve(Problem(Minimize(lambda_max(a))), solver = "ECOS"))
  expect_error(solve(Problem(Minimize(a)), solver = "SCS"))
})

test_that("Tests that a solver error is raised when a solver fails", {
  skip_on_cran()
  A <- matrix(rnorm(40*40), nrow = 40, ncol = 40)
  b <- A %*% matrix(rnorm(40))
  prob <- Problem(Minimize(sum_squares(A %*% Variable(40) - b)))
  expect_error(solve(prob, solver = "OSQP", max_iter = 1))
})

test_that("Tests problems with reshape_expr", {
  skip_on_cran()
  
  # Test on scalars.
  expect_equal(value(reshape_expr(1, c(1,1))), matrix(1))

  # Test vector to matrix.
  x <- Variable(4)
  mat <- cbind(c(1,-1), c(2,-2))
  vec <- matrix(1:4, ncol = 1)
  vec_mat <- cbind(c(1,2), c(3,4))
  expr <- reshape_expr(x, c(2,2))
  obj <- Minimize(sum_entries(mat %*% expr))
  prob <- Problem(obj, list(x == vec))
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, sum(mat %*% vec_mat), tolerance = TOL)

  # Test on matrix to vector.
  c <- 1:4
  expr <- reshape_expr(A, c(4,1))
  obj <- Minimize(t(expr) %*% c)
  constraints <- list(A == rbind(c(-1,-2), c(3,4)))
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, 20, tolerance = TOL)
  expect_equal(result$getValue(expr), matrix(c(-1,-2,3,4)))
  expect_equal(result$getValue(reshape_expr(expr, c(2,2))), rbind(c(-1,-2), c(3,4)))

  # Test matrix to matrix.
  expr <- reshape_expr(C, c(2,3))
  mat <- rbind(c(1,-1), c(2,-2))
  C_mat <- rbind(c(1,4), c(2,5), c(3,6))
  obj <- Minimize(sum_entries(mat %*% expr))
  prob <- Problem(obj, list(C == C_mat))
  result <- solve(prob, solver = "SCS")
  reshaped = matrix(C_mat, nrow = 2, ncol = 3, byrow = FALSE)
  expect_equal(result$value, sum(mat %*% reshaped), tolerance = TOL)
  expect_equal(result$getValue(expr), C_mat, tolerance = TOL)

  # Test promoted expressions.
  c <- cbind(c(1,-1), c(2,-2))
  expr <- reshape_expr(c * a, c(1,4))
  obj <- Minimize(expr %*% matrix(1:4))
  prob <- Problem(obj, list(a == 2))
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, -6, tolerance = TOL)
  expect_equal(result$getValue(expr), 2*matrix(c, nrow = 1), tolerance = TOL)

  expr <- reshape_expr(c * a, c(4,1))
  obj <- Minimize(t(expr) %*% matrix(1:4))
  prob <- Problem(obj, list(a == 2))
  result <- solve(prob, solver =)
  expect_equal(result$value, -6, tolerance = TOL)
  expect_equal(result$getValue(expr), 2*matrix(c, ncol = 1), tolerance = TOL)
})

test_that("test problems with cumsum", {
  tt <- Variable(5)
  prob <- Problem(Minimize(sum(tt)), list(cumsum(tt, 2) >= -0.0001))
  result <- solve(prob, solver = "SCS", eps = 1e-8)
  expect_equal(result$value, -0.0001, tolerance = TOL)
})

test_that("test problems with cummax", {
  tt <- Variable(5)
  prob <- Problem(Maximize(sum(tt)), list(cummax(tt, 2) <= c(1,2,3,4,5)))
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 15, tolerance = TOL)
})

test_that("Tests problems with vec", {
  skip_on_cran()
  c <- 1:4
  expr <- vec(A)
  obj <- Minimize(t(expr) %*% c)
  constraints <- list(A == rbind(c(-1,-2), c(3,4)))
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, 20, tolerance = TOL)
  expect_equal(result$getValue(expr), matrix(c(-1,-2,3,4)))
})

test_that("Test a problem with diag", {
  skip_on_cran()
  C <- Variable(3,3)
  obj <- Maximize(C[1,3])
  constraints <- list(diag(C) == 1, C[1,2] == 0.6, C[2,3] == -0.3, C == Variable(3, 3, PSD = TRUE))
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, 0.583151, tolerance = 1e-2)
})

test_that("Test presolve with parameters", {
  skip_on_cran()
  
  # Test with parameters.
  gamma <- Parameter(nonneg = TRUE)
  x <- Variable()
  obj <- Minimize(x)
  prob <- Problem(obj, list(gamma == 1, x >= 0))
  value(gamma) <- 0
  result <- solve(prob, solver = "SCS")
  expect_equal(result$status, INFEASIBLE)

  value(gamma) <- 1
  result <- solve(prob, solver = "SCS")
  expect_equal(result$status, OPTIMAL)
})

test_that("Test that expressions with parameters are updated properly", {
  x <- Variable()
  y <- Variable()
  x0 <- Parameter()
  xSquared <- x0*x0 + 2*x0*(x - x0)

  # Initial guess for x.
  value(x0) <- 2

  # Make the constraints x^2 - y == 0.
  g <- xSquared - y

  # Set up the problem.
  obj <- abs(x - 1)
  prob <- Problem(Minimize(obj), list(g == 0))
  expect_false(is_dpp(prob))
  expect_warning(solve(prob, solver = "SCS"))

  value(x0) <- 1
  expect_warning(result <- solve(prob, solver = "SCS"))
  expect_equal(result$getValue(g), 0, tolerance = TOL)

  # Test multiplication.
  prob <- Problem(Minimize(x0*x), list(x == 1))
  value(x0) <- 2
  expect_warning(solve(prob, solver = "SCS"))
  value(x0)
  expect_warning(result <- solve(prob, solver = "SCS"))
  expect_equal(result$value, 1, tolerance = 1e-2)
})

test_that("Test positive definite constraints", {
  skip_on_cran()
  C <- Variable(3,3)
  obj <- Maximize(C[1,3])
  constraints <- list(diag(C) == 1, C[1,2] == 0.6, C[2,3] == -0.3, C == t(C), C %>>% 0)
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, 0.583151, tolerance = 1e-2)

  C <- Variable(2,2)
  obj <- Maximize(C[1,2])
  constraints <- list(C == 1, C %>>% rbind(c(2,0), c(0,2)))
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$status, INFEASIBLE)

  C <- Variable(2, 2, symmetric = TRUE)
  obj <- Minimize(C[1,1])
  constraints <- list(C %<<% rbind(c(2,0), c(0,2)))
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$status, UNBOUNDED)
})

test_that("Test the duals of PSD constraints", {
  skip_on_cran()
  if("CVXOPT" %in% installed_solvers()) {
    # Test the dual values with CVXOPT.
    C <- Variable(2, 2, symmetric = TRUE, name = "C")
    obj <- Maximize(C[1,1])
    constraints <- list(C %<<% rbind(c(2,0), c(0,2)))
    prob <- Problem(obj, constraints)
    result <- solve(prob, solver = "CVXOPT")
    expect_equal(result$value, 2, tolerance = TOL)
    
    psd_constr_dual <- result$getDualValue(constraints[[1]])
    C <- Variable(2, 2, symmetric = TRUE, name = "C")
    X <- Variable(2, 2, PSD = TRUE)
    obj <- Maximize(C[1,1])
    constraints <- list(X == rbind(c(2,0), c(0,2)) - C)
    prob <- Problem(obj, constraints)
    result <- solve(prob, solver = "CVXOPT")
    # Symmetrizing is valid, because the dual variable is with respect to an
    # unstructured equality constraint. Dual optimal solutions are non-unique
    # in this formulation, up to symmetrizing the dual variable.
    new_constr_dual <- (result$getDualValue(constraints[[1]]) + result$getDualValue(constraints[[2]]))/2
    expect_equal(new_constr_dual, psd_constr_dual, tolerance = TOL)
  }
  
  # Test dual values with SCS.
  C <- Variable(2, 2, symmetric = TRUE)
  obj <- Maximize(C[1,1])
  constraints <- list(C %<<% rbind(c(2,0), c(0,2)))
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, 2, tolerance = 1e-4)

  psd_constr_dual <- result$getDualValue(constraints[[1]])
  C <- Variable(2, 2, symmetric = TRUE)
  X <- Variable(2, 2, PSD = TRUE)
  obj <- Maximize(C[1,1])
  constraints <- list(X == rbind(c(2,0), c(0,2)) - C)
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$getDualValue(constraints[[1]]), psd_constr_dual, tolerance = TOL)

  # Test dual values with SCS that have off-diagonal entries.
  C <- Variable(2, 2, symmetric = TRUE)
  obj <- Maximize(C[1,2] + C[2,1])
  constraints <- list(C %<<% rbind(c(2,0), c(0,2)), C >= 0)
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, 4, tolerance = 1e-3)

  psd_constr_dual <- result$getDualValue(constraints[[1]])
  C <- Variable(2, 2, symmetric = TRUE)
  X <- Variable(2, 2, PSD = TRUE)
  obj <- Maximize(C[1,2] + C[2,1])
  constraints <- list(X == rbind(c(2,0), c(0,2)) - C, C >= 0)
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$getDualValue(constraints[[1]]), psd_constr_dual, tolerance = 1e-3)
})

test_that("Test geo_mean function", {
  skip_on_cran()
  x <- Variable(2)
  cost <- geo_mean(x)
  prob <- Problem(Maximize(cost), list(x <= 1))
  result <- solve(prob, solver = "SCS", eps = 1e-5)
  expect_equal(result$value, 1, tolerance = TOL)

  prob <- Problem(Maximize(cost), list(sum(x) <= 1))
  result <- solve(prob, solver = "SCS", eps = 1e-5)
  expect_equal(result$getValue(x), matrix(c(0.5,0.5)), tolerance = TOL)

  x <- Variable(3,3)
  expect_error(geo_mean(x))

  x <- Variable(3,1)
  g <- geo_mean(x)
  # expect_equal(g@w, gmp::as.bigq(1,3)*3, tolerance = TOL)

  x <- Variable(1,5)
  g <- geo_mean(x)
  # expect_equal(g@w, gmp::as.bigq(1,5)*5, tolerance = TOL)

  # Check that we get the right answer for max geo_mean(x) s.t. sum(x) <= 1.
  p <- c(0.07, 0.12, 0.23, 0.19, 0.39)

  short_geo_mean <- function(x, p) {
      p <- as.numeric(p)/sum(p)
      x <- as.numeric(x)
      return(prod(x^p))
  }

  x <- Variable(5)
  prob <- Problem(Maximize(geo_mean(x, p)), list(sum(x) <= 1))
  result <- solve(prob, solver = "SCS", eps = 1e-5)
  x <- as.numeric(result$getValue(x))
  x_true <- p/sum(p)

  expect_true(is.allclose(result$value, result$getValue(geo_mean(x, p))))
  expect_true(is.allclose(result$value, short_geo_mean(x, p)))
  expect_true(is.allclose(x, x_true, 1e-3))
  
  # Check that we get the right answer for max geo_mean(x) s.t. p_norm(x) <= 1.
  x <- Variable(5)
  prob <- Problem(Maximize(geo_mean(x, p)), list(p_norm(x) <= 1))
  result <- solve(prob, solver = "SCS", eps = 1e-5)
  x <- as.numeric(result$getValue(x))
  x_true <- sqrt(p/sum(p))

  expect_true(is.allclose(result$value, result$getValue(geo_mean(x, p))))
  expect_true(is.allclose(result$value, short_geo_mean(x, p)))
  expect_true(is.allclose(x, x_true, 1e-3))

  # The following 3 tests check vstack and hstack input to geo_mean
  # The following 3 formulations should be equivalent
  n <- 5
  x_true <- matrix(rep(1,n))
  x <- Variable(n)

  result <- solve(Problem(Maximize(geo_mean(x)), list(x <= 1)), solver = "SCS")
  xval <- as.numeric(result$getValue(x))
  expect_true(is.allclose(xval, x_true, 1e-3))

  args <- lapply(seq_len(n), function(ind) x[ind])
  y <- do.call(vstack, args)
  result <- solve(Problem(Maximize(geo_mean(y)), list(x <= 1)), solver = "SCS")
  xval <- as.numeric(result$getValue(x))
  expect_true(is.allclose(xval, x_true, 1e-3))

  y <- do.call(hstack, args)
  result <- solve(Problem(Maximize(geo_mean(y)), list(x <= 1)), solver = "SCS")
  xval <- as.numeric(result$getValue(x))
  expect_true(is.allclose(xval, x_true, 1e-3))
})

test_that("Test p_norm function", {
  skip_on_cran()
  x <- Variable(3, name = "x")
  avec <- c(1.0, 2, 3)

  # TODO: Add -1, 0.5, 0.3, -2.3 and testing positivity constraints

  for(p in c(1, 1.6, 1.3, 2, 1.99, 3, 3.7, Inf)) {
    prob <- Problem(Minimize(p_norm(x, p = p)), list(t(x) %*% avec >= 1))
    result <- solve(prob, solver = "ECOS", verbose = TRUE)

    # Formula is true for any a >= 0 with p > 1
    if(p == Inf)
      x_true <- rep(1, length(avec))/sum(avec)
    else if(p == 1) {
      # Only works for the particular a = c(1,2,3)
      x_true <- c(0,0,1.0/3)
    } else
      x_true <- avec^(1.0/(p-1)) / as.numeric(avec %*% (avec^(1.0/(p-1))))

    x_alg <- as.vector(result$getValue(x))
    print(paste("p =", p))
    expect_true(is.allclose(x_alg, x_true, 1e-2))
    expect_true(is.allclose(result$value, .p_norm(x_alg, p)))
    expect_true(is.allclose(.p_norm(x_alg, p), result$getValue(p_norm(x_alg, p))))
  }
})

test_that("Test p_norm concave", {
  skip_on_cran()
  x <- Variable(3, name = "x")

  # Test positivity constraints
  a <- c(-1.0, 2, 3)
  for(p in c(-1, 0.5, 0.3, -2.3)) {
    prob <- Problem(Minimize(sum_entries(abs(x-a))), list(p_norm(x,p) >= 0))
    result <- solve(prob, solver = "ECOS")
    expect_true(is.allclose(result$value, 1))
  }

  a <- c(1.0, 2, 3)
  for(p in c(-1, 0.5, 0.3, -2.3)) {
    prob <- Problem(Minimize(sum_entries(abs(x-a))), list(p_norm(x,p) >= 0))
    result <- solve(prob, solver = "ECOS")
    expect_equal(result$value, 0, tolerance = 1e-6)
  }
})

test_that("Test power function", {
  skip_on_cran()
  x <- Variable()
  prob <- Problem(Minimize(power(x, 1.7) + power(x, -2.3) - power(x, 0.45)))
  result <- solve(prob, solver = "SCS", eps = 1e-5)
  xs <- result$getValue(x)
  expect_true(abs(1.7*xs^0.7 - 2.3*xs^-3.3 - 0.45*xs^-0.55) <= 1e-3)
})

test_that("Test a problem with multiply by a scalar", {
  skip_on_cran()
  Tnum <- 10
  Jnum <- 20
  rvec <- matrix(stats::rnorm(Tnum*Jnum), nrow = Tnum, ncol = Jnum)
  dy <- matrix(stats::rnorm(2*Tnum), nrow = 2*Tnum, ncol = 1)
  theta <- Variable(Jnum)

  delta <- 1e-3
  loglambda <- rvec %*% theta   # rvec: TxJ regressor matrix, theta: (Jx1) cvx variable
  a <- multiply(dy[1:Tnum], loglambda)  # size (Tx1)
  b1 <- exp(loglambda)
  b2 <- multiply(delta, b1)
  cost <- -a + b1

  cost <- -a + b2  # size (Tx1)
  prob <- Problem(Minimize(sum_entries(cost)))
  result <- solve(prob, solver = "SCS")

  obj <- Minimize(sum_entries(multiply(2, x)))
  prob <- Problem(obj, list(x == 2))
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, 8, tolerance = TOL)
})

test_that("Test bug with 64-bit integers", {
  skip_on_cran()
  if(!require(bit64)) {
    print("bit64 library not found. Skipping test.")
    return()
  }
  
  q <- Variable(as.integer64(2))
  objective <- Minimize(norm1(q))
  problem <- Problem(objective)
  result <- solve(problem, solver = "SCS")
  print(result$getValue(q))
})

test_that("Test bug with negative slice", {
  skip_on_cran()
  x <- Variable(2)
  objective <- Minimize(x[1] + x[2])
  constraints <- list(x[-2:nrow(x)] >= 1)
  problem <- Problem(objective, constraints)
  result <- solve(problem, solver = "SCS", eps = 1e-6)
  expect_equal(result$getValue(x), matrix(c(1,1)), tolerance = TOL)
})

test_that("Test p_norm with axis != 2", {
  skip_on_cran()
  b <- 0:1
  X <- Variable(2, 10)
  expr <- p_norm(X, p = 2, axis = 1) - b
  con <- list(expr <= 0)
  obj <- Maximize(sum(X))
  prob <- Problem(obj, con)
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$getValue(expr), matrix(c(0,0)), tolerance = TOL)
  
  b <- 0:9
  X <- Variable(10, 2)
  expr <- p_norm(X, p = 2, axis = 1) - b
  con <- list(expr <= 0)
  obj <- Maximize(sum(X))
  prob <- Problem(obj, con)
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$getValue(expr), matrix(rep(0, 10)), tolerance = TOL)
})

test_that("Test constraints that evaluate to booleans", {
  skip_on_cran()
  x <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x), list(TRUE))
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$getValue(x), 0, tolerance = TOL)
  
  prob <- Problem(Minimize(x), as.list(rep(TRUE, 10)))
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$getValue(x), 0, tolerance = TOL)
  
  prob <- Problem(Minimize(x), c(list(42 <= x), as.list(rep(TRUE, 10))))
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$getValue(x), 42, tolerance = TOL)
  
  prob <- Problem(Minimize(x), c(list(TRUE), list(42 <= x), as.list(rep(TRUE, 10))))
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$getValue(x), 42, tolerance = TOL)
  
  prob <- Problem(Minimize(x), list(FALSE))
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$status, INFEASIBLE)
  
  prob <- Problem(Minimize(x), as.list(rep(FALSE, 10)))
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$status, INFEASIBLE)
  
  prob <- Problem(Minimize(x), c(as.list(rep(TRUE, 10)), list(FALSE), as.list(rep(TRUE, 10))))
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$status, INFEASIBLE)
  
  # Only TRUEs, but infeasible solution since x must be non-negative.
  prob <- Problem(Minimize(x), c(list(TRUE), list(x <= -42), as.list(rep(TRUE, 10))))
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$status, INFEASIBLE)
})

test_that("Test a problem with an invalid constraint", {
  skip_on_cran()
  x <- Variable()
  expect_error(Problem(Minimize(x), list(sum(x))), "Problem has an invalid constraint")
})

test_that("Test the pos and neg attributes", {
  skip_on_cran()
  x <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x))
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$getValue(x), 0, tolerance = TOL)
  
  x <- Variable(neg = TRUE)
  prob <- Problem(Maximize(x))
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$getValue(x), 0, tolerance = TOL)
})

test_that("Test saving and loading problems", {
  skip_on_cran()
  prob <- Problem(Minimize(2*a + 3), list(a >= 1))
  saveRDS(prob, file = "prob_test.RDS")   # TODO: Figure out where to save this in unit test runs (esp. on CRAN).
  new_prob <- loadRDS("prob_test.RDS")
  result <- solve(new_prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 5.0, tolerance = TOL)
  expect_equal(result$getValue(variables(new_prob)[[1]]), 1.0, tolerance = TOL)
})

test_that("Test problem with sparse int8 matrix", {
  skip_on_cran()
  a <- Variable(3, 1)
  q <- matrix(c(1.88922129, 0.06938685, 0.91948919), ncol = 1)
  P <- rbind(c(280.64, -49.84, -80),
             c(-49.84, 196.04, 139),
             c(-80, 139, 106))
  D_dense <- rbind(c(-1, 1, 0,  0, 0, 0),
                   c(0, -1, 1,  0, 0, 0),
                   c(0,  0, 0, -1, 1, 0))
  D_dense <- as.integer(D_dense)
  D_sparse <- Matrix(D_dense, sparse = TRUE)
  
  make_problem <- function(D) {
    obj <- Minimize(0.5*quad_form(a, P), - t(a) %*% q)
    expect_true(is_dcp(obj))
    
    alpha <- Parameter(nonneg = TRUE, value = 2)
    constraints <- list(a >= 0, -alpha <= t(D) %*% a, t(D) %*% a <= alpha)
    
    prob <- Problem(obj, constraints)
    result <- solve(prob, solver = "ECOS")
    expect_equal(result$status, "optimal")
    return(list(prob, result))
  }
  
  expected_coef <- matrix(c(-0.011728003147, 0.011728002895, 0.000000000252,
                            -0.017524801335, 0.017524801335, 0.), nrow = 1)
  
  tmp <- make_problem(D_dense)
  prob <- tmp[[1]]
  result <- tmp[[2]]
  coef_dense <- t(result$getValue(a)) %*% D_dense
  expect_equal(expected_coef, coef_dense, tolerance = TOL)
  
  tmp <- make_problem(D_sparse)
  prob <- tmp[[1]]
  result <- tmp[[2]]
  coef_sparse <- t(result$getValue(a)) %*% D_sparse
  expect_equal(expected_coef, coef_sparse, tolerance = TOL)
})

test_that("Test QP code path with special indexing", {
  skip_on_cran()
  
  # TODO: Is this test necessary given how R handles indexing? (No slice object like Python).
  x <- Variable(1, 3)
  y <- sum(x[,1:2], axis = 1)
  cost <- QuadForm(y, diag(1))
  prob <- Problem(Minimize(cost))
  result1 <- solve(prob, solver = "SCS")
  
  x <- Variable(1, 3)
  y <- sum(x[,c(1,2)], axis = 1)
  cost <- QuadForm(y, diag(1))
  prob <- Problem(Minimize(cost))
  result2 <- solve(prob, solver = "SCS")
  expect_equal(result1$value, result2$value, tolerance = TOL)
})

test_that("Test a problem with indicators", {
  skip_on_cran()
  n <- 5
  m <- 2
  q <- 0:(n-1)
  a <- matrix(1, nrow = m, ncol = n)
  b <- matrix(1, nrow = m, ncol = 1)
  x <- Variable(n, 1, name = "x")
  
  constraints <- list(a %*% x == b)
  objective <- Minimize((1/2)*square(t(q) %*% x) + indicator(constraints))
  problem <- Problem(objective)
  result1 <- solve(problem, solver = "SCS", eps = 1e-5)
  solution1 <- result1$value
  
  # Without indicators.
  objective <- Minimize((1/2)*square(t(q) %*% x))
  problem <- Problem(objective, constraints)
  result2 <- solve(problem, solver = "SCS", eps = 1e-5)
  solution2 <- result2$value
  expect_equal(solution1, solution2, tolerance = TOL)
})

test_that("Test that rmul works with 1x1 matrices", {
  skip_on_cran()
  x <- matrix(c(4144.30127531))
  y <- matrix(c(7202.52114311))
  z <- Variable(1, 1)
  objective <- Minimize(quad_form(z, x) - 2*t(z) %*% y)
  
  prob <- Problem(objective)
  res <- solve(prob, solver = "OSQP", verbose = TRUE)
  result1 <- res$value
  
  x <- 4144.30127531
  y <- 7202.52114311
  z <- Variable()
  objective <- Minimize(x*z^2 - 2*z*y)
  
  prob <- Problem(objective)
  res <- solve(prob, solver = "OSQP", verbose = TRUE)
  expect_equal(res$value, result1, tolerance = TOL)
})

test_that("Test reshape of a min with axis = 2", {
  skip_on_cran()
  x <- Variable(5, 2)
  y <- Variable(5, 2)
  
  stacked_flattened <- vstack(vec(x), vec(y))   # (2, 10).
  minimum <- min_entries(stacked_flatten, axis = 2)   # (10,)
  reshaped_minimum <- reshape_expr(minimum, c(5, 2))   # (5, 2).
  
  obj <- sum(reshaped_minimum)
  problem <- Problem(Maximize(obj), list(x == 1, y == 2))
  result <- solve(problem, solver = "SCS")
  expect_equal(result$value, 10, tolerance = TOL)
})

test_that("Test a problem with constant values only that is infeasible", {
  skip_on_cran()
  p <- Problem(Maximize(0), list(Constant(0) == 1))
  result <- solve(p, solver = "SCS")
  expect_equal(result$status, INFEASIBLE)
})

test_that("Test Huber regression works with SCS", {
  skip_on_cran()
  
  # See CVXPY issue #1370.
  set.seed(1)
  m <- 5
  n <- 2
  
  x0 <- matrix(rnorm(n), nrow = n, ncol = 1)
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  b <- A %*% x0 + 0.01*matrix(rnorm(m))
  
  # Add outlier noise.
  k <- as.integer(0.02*m)
  idx <- sample.int(m, size = k)
  b[idx] <- b[idx] + 10*matrix(rnorm(k))
  
  x <- Variable(n)
  prob <- Problem(Minimize(sum(huber(A %*% x - b))))
  result <- solve(prob, solver = "SCS")
})

test_that("Test a complex rmul expression with a parameter", {
  skip_on_cran()
  
  # See CVXPY issue #1555.
  b <- Variable(1)
  param <- Parameter(1)
  
  constraints <- list()
  objective <- Minimize((2*b) %*% param)
  prob <- Problem(objective, constraints)
  
  value(param) <- matrix(1)
  result <- solve(prob)
  expect_equal(result$value, -Inf)
})

test_that("Test the cumsum axis bug with row or column matrix", {
  skip_on_cran()
  
  # See CVXPY issue #1678.
  n <- 5
  
  # Solve for axis = 1.
  x1 <- Variable(n, 1)
  expr1 <- cumsum(x1, axis = 1)
  prob1 <- Problem(Minimize(0), list(expr1 == 1))
  result <- solve(prob1)
  expect <- matrix(1, nrow = n, ncol = 1)
  expect_equal(result$getValue(expr1), expect, tolerance = TOL)
  
  # Solve for axis = 2.
  x2 <- Variable(1, n)
  expr2 <- cumsum(x2, axis = 2)
  prob2 <- Problem(Minimize(0), list(expr2 == 1))
  result <- solve(prob2)
  expect <- matrix(1, nrow = 1, ncol = n)
  expect_equal(result$getValue(expr2), expect, tolerance = TOL)
})

test_that("Test the cummax axis bug with row or column matrix", {
  skip_on_cran()
  
  # See CVXPY issue #1678.
  n <- 5
  
  # Solve for axis = 1.
  x1 <- Variable(n, 1)
  expr1 <- cummax(x1, axis = 1)
  prob1 <- Problem(Maximize(sum(x1)), list(expr1 == 1))
  result <- solve(prob1)
  expect <- matrix(1, nrow = n, ncol = 1)
  expect_equal(result$getValue(expr1), expect, tolerance = TOL)
  
  # Solve for axis = 2.
  x2 <- Variable(1, n)
  expr2 <- cummax(x2, axis = 2)
  prob2 <- Problem(Maximize(sum(x2)), list(expr2 == 1))
  result <- solve(prob2)
  expect <- matrix(1, nrow = 1, ncol = n)
  expect_equal(result$getValue(expr2), expect, tolerance = TOL)
})
