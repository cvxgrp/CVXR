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

test_that("test the variables method", {
  skip_on_cran()
  p <- Problem(Minimize(a), list(a <= x, b <= A + 2))
  vars_ <- variables(p)
  ref <- list(a, x, b, A)
  mapply(function(v, r) { expect_equal(v, r) }, vars_, ref)
})

test_that("test the parameters method", {
  skip_on_cran()
  p1 <- Parameter()
  p2 <- Parameter(3, nonpos = TRUE)
  p3 <- Parameter(4, 4, nonneg = TRUE)
  p <- Problem(Minimize(p1), list(a + p1 <= p2, b <= p3 + p3 + 2))
  params <- parameters(p)
  ref <- c(p1, p2, p3)
  mapply(function(p, r) { expect_equal(p, r) }, params, ref)
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

  # Single scalar constants
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
  ref <- size(p1) + size(p2) + size(p3) + length(constants)  # 2 and c2 %*% c1 are both single scalar constants
  expect_equal(n_data, ref)

  # num_scalar_eq_constr
  n_eq_constr <- size_metrics(p)@num_scalar_eq_constr
  ref <- prod(dim(c2 %*% c1))
  expect_equal(n_eq_constr, ref)

  # num_scalar_leq_constr
  n_leq_constr <- size_metrics(p)@num_scalar_leq_constr
  ref <- size(p3) + size(p2)
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
  expect_true(result$solve_time > 0)
  expect_true(result$setup_time > 0)
  expect_true(result$num_iters > 0)
})

test_that("Test the get_problem_data method", {
  skip_on_cran()
  data <- get_problem_data(Problem(Minimize(exp(a) + 2)), "SCS")[[1]]
  dims <- data[[ConicSolver()@dims]]
  expect_equal(dims@exp, 1)
  expect_equal(length(data[["c"]]), 2)
  expect_equal(dim(data[["A"]]), c(3,2))

  data <- get_problem_data(Problem(Minimize(p_norm(x) + 3)), "ECOS")[[1]]
  dims <- data[[ConicSolver()@dims]]
  expect_equal(dims@soc, list(3))
  expect_equal(length(data[["c"]]), 3)
  expect_equal(dim(data[["A"]]), c(0,0))
  expect_equal(dim(data[["G"]]), c(3,3))

  if("CVXOPT" %in% INSTALLED_SOLVERS) {
    data <- get_problem_data(Problem(Minimize(p_norm(x) + 3)), "CVXOPT")[[1]]
    dims <- data[["dims"]]
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
  data <- list(c = args[["c"]], A = args[["A"]], b = args[["b"]])
  cones <- SCS.dims_to_solver_dict(args[[ConicSolver()@dims]])
  solution <- scs::scs(data$A, data$b, data$c, cones)
  prob <- Problem(Minimize(exp(a)), list(a == 0))
  result <- unpack_results(prob, solution, chain, inv)
  expect_equal(result$getValue(a), 0, tolerance = 1e-3)
  expect_equal(result$value, 1, tolerance = 1e-3)
  expect_equal(result$status, "optimal")

  prob <- Problem(Minimize(p_norm(x)), list(x == 0))
  tmp <- get_problem_data(prob, solver = "ECOS")
  args <- tmp[[1]]
  chain <- tmp[[2]]
  inv <- tmp[[3]]
  cones <- ECOS.dims_to_solver_dict(args[[ConicSolver()@dims]])
  solution <- ECOSolveR::ECOS_csolve(args[["c"]], args[["G"]], args[["h"]], cones, args[["A"]], args[["b"]])
  prob <- Problem(Minimize(p_norm(x)), list(x == 0))
  result <- unpack_results(prob, solution, chain, inv)
  expect_equal(result$getValue(x), matrix(c(0,0)))
  expect_equal(result$value, 0)
  expect_equal(result$status, "optimal")
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
  # Vector variables
  p <- Problem(Maximize(x[1,1]), list(x[1,1] <= 2, x[2,1] == 3))
  result <- solve(p)
  expect_equal(result$value, 2, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(2,3)))

  n <- 10
  Aloc <- matrix(0:(n^2-1), nrow = n, ncol = n)
  xloc <- Variable(n,n)
  p <- Problem(Minimize(sum_entries(xloc)), list(xloc == Aloc))
  result <- solve(p)
  answer <- n*n*(n*n+1)/2 - n*n
  expect_equal(result$value, answer)

  # Matrix variables
  obj <- A[1,1] + A[1,2] + A[2,2] + A[2,1]
  p <- Problem(Maximize(obj), list(A <= rbind(c(1,-2), c(-3,4))))
  result <- solve(p)
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(A), rbind(c(1,-2), c(-3,4)))

  # Indexing arithmetic expressions
  exp <- cbind(c(1,2), c(3,4)) %*% z + x
  p <- Problem(Minimize(exp[2,1]), list(x == z, z == c(1,2)))
  result <- solve(p)
  expect_equal(result$value, 12, tolerance = TOL)
  expect_equal(result$getValue(x), result$getValue(z), tolerance = TOL)
})

test_that("Test problems with slicing", {
  skip_on_cran()
  p <- Problem(Maximize(sum_entries(C)), list(C[2:3,] <= 2, C[1,] == 1))
  result <- solve(p)
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equal(result$getValue(C), cbind(c(1,2,2), c(1,2,2)))

  p <- Problem(Maximize(sum_entries(C[seq(1,3,2),2])), list(C[2:3,] <= 2, C[1,] == 1))
  result <- solve(p)
  expect_equal(result$value, 3, tolerance = TOL)
  expect_equal(result$getValue(C[seq(1,3,2),2]), matrix(c(1,2)))

  p <- Problem(Maximize(sum_entries((C[1:2,] + A)[,1:2])), list(C[2:3,] <= 2, C[1,] == 1, (A + B)[,1] == 3,
                                                             (A + B)[,2] == 2, B == 1))
  result <- solve(p)
  expect_equal(result$value, 12, tolerance = TOL)
  expect_equal(result$getValue(C[1:2,]), cbind(c(1,2), c(1,2)), tolerance = TOL)
  expect_equal(result$getValue(A), cbind(c(2,2),c(1,1)), tolerance = TOL)

  p <- Problem(Maximize(matrix(c(3,4), nrow = 1, ncol = 2) %*% (C[1:2,] + A)[,1]),
               list(C[2:3,] <= 2, C[1,] == 1, matrix(c(1,2), nrow = 1, ncol = 2) %*% (A + B)[,1] == 3,
                    (A + B)[,2] == 2, B == 1, 3*A[,1] <= 3))
  result <- solve(p)
  expect_equal(result$value, 12, tolerance = TOL)
  expect_equal(result$getValue(C[1:2,1]), matrix(c(1,2)), tolerance = TOL)
  expect_equal(result$getValue(A), cbind(c(1,-0.5), c(1,1)), tolerance = TOL)

  p <- Problem(Minimize(norm2((C[1:2,] + A)[,1])), list(C[2:3,] <= 2, C[1,] == 1, (A + B)[,1] == 3, (A + B)[,2] == 2, B == 1))
  result <- solve(p)
  expect_equal(result$value, 3, tolerance = TOL)
  expect_equal(result$getValue(C[1:2,1]), matrix(c(1,-2)), tolerance = 1e-4)
  expect_equal(result$getValue(A), cbind(c(2,2), c(1,1)), tolerance = TOL)

  # Transpose of slice
  p <- Problem(Maximize(sum_entries(C)), list(t(C[2:3,]) <= 2, t(C[1,]) == 1))
  result <- solve(p)
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equal(result$getValue(C), cbind(c(1,2,2), c(1,2,2)))
})

test_that("Test the vstack function", {
  skip_on_cran()
  c <- matrix(1, nrow = 1, ncol = 5)
  p <- Problem(Minimize(c %*% vstack(x, y)), list(x == c(1,2), y == c(3,4,5)))
  result <- solve(p)
  expect_equal(result$value, 15, tolerance = TOL)

  c <- matrix(1, nrow = 1, ncol = 4)
  p <- Problem(Minimize(c %*% vstack(x, x)), list(x == c(1,2)))
  result <- solve(p)
  expect_equal(result$value, 6, tolerance = TOL)

  c <- matrix(1, nrow = 2, ncol = 2)
  p <- Problem(Minimize(sum_entries(vstack(A, C))), list(A >= 2*c, C == -2))
  result <- solve(p)
  expect_equal(result$value, -4, tolerance = TOL)

  c <- matrix(1, nrow = 1, ncol = 2)
  p <- Problem(Minimize(sum_entries(vstack(c %*% A, c %*% B))), list(A >= 2, B == -2))
  result <- solve(p)
  expect_equal(result$value, 0, tolerance = TOL)

  c <- matrix(c(1,-1), nrow = 2, ncol = 1)
  p <- Problem(Minimize(t(c) %*% vstack(a^2, sqrt(b))), list(a == 2, b == 16))
  expect_error(solve(p))
})

test_that("Test the hstack function", {
  skip_on_cran()
  c <- matrix(1, nrow = 1, ncol = 5)
  p <- Problem(Minimize(c %*% t(hstack(t(x), t(y)))), list(x == c(1,2), y == c(3,4,5)))
  result <- solve(p)
  expect_equal(result$value, 15, tolerance = TOL)

  c <- matrix(1, nrow = 1, ncol = 4)
  p <- Problem(Minimize(c %*% t(hstack(t(x), t(x)))), list(x == c(1,2)))
  result <- solve(p)
  expect_equal(result$value, 6, tolerance = TOL)

  c <- matrix(1, nrow = 2, ncol = 2)
  p <- Problem(Minimize(sum_entries(hstack(t(A), t(C)))), list(A >= 2*c, C == -2))
  result <- solve(p)
  expect_equal(result$value, -4, tolerance = TOL)

  D <- Variable(3,3)
  expr <- hstack(C, D)
  p <- Problem(Minimize(expr[1,2] + sum_entries(hstack(expr, expr))), list(C >= 0, D >= 0, D[1,1] == 2, C[1,2] == 3))
  result <- solve(p)
  expect_equal(result$value, 13, tolerance = TOL)

  c <- matrix(c(1,-1), nrow = 2, ncol = 1)
  p <- Problem(Minimize(t(c) %*% t(hstack(t(a^2), t(sqrt(b))))), list(a == 2, b == 16))
  expect_error(solve(p))
})

test_that("Test using a CVXR expression as an objective", {
  skip_on_cran()
  expect_error(Problem(x+2))
})

test_that("Test variable transpose", {
  skip_on_cran()
  p <- Problem(Minimize(sum_entries(x)), list(t(x) >= matrix(c(1,2), nrow = 1, ncol = 2)))
  result <- solve(p)
  expect_equal(result$value, 3, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(1,2)))

  p <- Problem(Minimize(sum_entries(C)), list(matrix(c(1,1), nrow = 1, ncol = 2) %*% t(C) >= matrix(0:2, nrow = 1, ncol = 3)))
  result <- solve(p)
  value <- result$getValue(C)

  constraints <- lapply(1:3, function(i) { 1*C[i,1] + 1*C[i,2] >= (i-1) })
  p <- Problem(Minimize(sum_entries(C)), constraints)
  result2 <- solve(p)
  expect_equal(result$value, result2$value)
  expect_equal(result$getValue(C), value)

  p <- Problem(Minimize(A[1,2] - t(A)[2,1]), list(A == cbind(c(1,2), c(3,4))))
  result <- solve(p)
  expect_equal(result$value, 0, tolerance = TOL)

  exp <- t(-x)
  p <- Problem(Minimize(sum_entries(x)), list(t(-x) <= 1))
  result <- solve(p)
  expect_equal(result$value, -2, tolerance = TOL)

  c <- matrix(c(1,-1), nrow = 2, ncol = 1)
  p <- Problem(Minimize(max_elemwise(t(c), 2, 2 + t(c))[2]))
  result <- solve(p)
  expect_equal(result$value, 2, tolerance = TOL)

  c <- cbind(c(1,-1,2), c(1,-1,2))
  p <- Problem(Minimize(sum_entries(t(max_elemwise(c, 2, 2+c))[,1])))
  result <- solve(p)
  expect_equal(result$value, 6, tolerance = TOL)

  c <- cbind(c(1,-1,2), c(1,-1,2))
  p <- Problem(Minimize(sum_entries(t(t(c)^2)[,1])))
  result <- solve(p)
  expect_equal(result$value, 6, tolerance = TOL)

  # Slice of transpose
  p <- Problem(Maximize(sum_entries(C)), list(t(C)[,2:3] <= 2, t(C)[,1] == 1))
  result <- solve(p)
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equal(result$getValue(C), cbind(c(1,2,2), c(1,2,2)))
})

test_that("Test multiplication on the left by a non-constant", {
  skip_on_cran()
  c <- matrix(c(1,2), nrow = 2, ncol = 1)
  p <- Problem(Minimize(t(c) %*% A %*% c), list(A >= 2))
  result <- solve(p)
  expect_equal(result$value, 18, tolerance = TOL)

  p <- Problem(Minimize(a*2), list(a >= 2))
  result <- solve(p)
  expect_equal(result$value, 4, tolerance = TOL)

  p <- Problem(Minimize(t(x) %*% c), list(x >= 2))
  result <- solve(p)
  expect_equal(result$value, 6, tolerance = TOL)

  p <- Problem(Minimize((t(x) + t(z)) %*% c), list(x >= 2, z >= 1))
  result <- solve(p)
  expect_equal(result$value, 9, tolerance = TOL)
})

test_that("Test redundant constraints", {
  skip_on_cran()
  obj <- Minimize(sum_entries(x))
  constraints <- list(x == 2, x == 2, t(x) == 2, x[1] == 2)
  p <- Problem(obj, constraints)
  result <- solve(p, solver = "ECOS")
  expect_equal(result$value, 4, tolerance = TOL)

  obj <- Minimize(sum_entries(x^2))
  constraints <- list(x == x)
  p <- Problem(obj, constraints)
  result <- solve(p, solver = "ECOS")
  expect_equal(result$value, 0, tolerance = TOL)
})

test_that("Test that symmetry is enforced", {
  skip_on_cran()
  p <- Problem(Minimize(lambda_max(A)), list(A >= 2))
  result <- solve(p)
  expect_equal(result$getValue(A), t(result$getValue(A)), tolerance = 1e-3)

  p <- Problem(Minimize(lambda_max(A)), list(A == cbind(c(1,2), c(3,4))))
  result <- solve(p)
  expect_equal(result$status, "infeasible")
})

test_that("Test SDP", {
  skip_on_cran()
  # Ensure SDP constraints enforce transpose
  obj <- Maximize(A[2,1] - A[1,2])
  p <- Problem(obj, list(lambda_max(A) <= 100, A[1,1] == 2, A[2,2] == 2, A[2,1] == 2))
  result <- solve(p)
  expect_equal(result$value, 0, tolerance = 1e-3)
})

test_that("Test getting values for expressions", {
  skip_on_cran()
  diff_exp <- x - z
  inf_exp <- norm_inf(diff_exp)
  sum_entries_exp <- 5 + p_norm(z,1) + p_norm(x,1) + inf_exp
  constr_exp <- norm2(x + z)
  obj <- norm2(sum_entries_exp)
  p <- Problem(Minimize(obj), list(x >= c(2,3), z <= c(-1,-4), constr_exp <= 2))
  result <- solve(p)
  expect_equal(result$value, 22, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(2,3)), tolerance = TOL)
  expect_equal(result$getValue(z), matrix(c(-1,-4)), tolerance = TOL)

  # Expression values
  xs <- result$getValue(x)
  zs <- result$getValue(z)
  expect_equal(result$getValue(diff_exp), xs - zs, tolerance = TOL)
  expect_equal(result$getValue(inf_exp), norm(xs - zs, "I"))
  expect_equal(result$getValue(sum_entries_exp), 5 + norm(zs, "1") + norm(xs, "1") + norm(xs - zs, "I"))
  expect_equal(result$getValue(constr_exp), base::norm(xs + zs, "2"))
  expect_equal(result$getValue(obj), result$value)
})

test_that("Test multiplication by zero", {
  skip_on_cran()
  exp <- 0*a
  expect_equal(value(exp), matrix(0))
  obj <- Minimize(exp)
  p <- Problem(obj)
  result <- solve(p)
  expect_equal(result$value, 0, tolerance = TOL)
  expect_false(is.na(result$getValue(a)))
})

test_that("Tests a problem with division", {
  skip_on_cran()
  obj <- Minimize(norm_inf(A/5))
  p <- Problem(obj, list(A >= 5))
  result <- solve(p)
  expect_equal(result$value, 1, tolerance = TOL)
})

test_that("Tests problems with multiply", {
  skip_on_cran()
  c <- cbind(c(1,-1), c(2,-2))
  expr <- multiply(c, A)
  obj <- Minimize(norm_inf(expr))
  p <- Problem(obj, list(A == 5))
  result <- solve(p)
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equal(result$getValue(expr), cbind(c(5,-5), c(10,-10)), tolerance = TOL)

  # Test with a sparse matrix
  c <- Matrix::sparseMatrix(i = c(1,2), j = c(1,1), x = c(1,2))
  expr <- multiply(c, x)
  obj <- Minimize(norm_inf(expr))
  p <- Problem(obj, list(x == 5))
  result <- solve(p)
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equivalent(as.matrix(result$getValue(expr)), matrix(c(5,10)))

  # Test promotion
  c <- cbind(c(1,-1), c(2,-2))
  expr <- multiply(c, a)
  obj <- Minimize(norm_inf(expr))
  p <- Problem(obj, list(a == 5))
  result <- solve(p)
  expect_equal(result$value, 10, tolerance = TOL)
  expect_equivalent(as.matrix(result$getValue(expr)), cbind(c(5,-5), c(10,-10)))
})

test_that("Tests that errors occur when you use an invalid solver", {
  skip_on_cran()
  expect_error(solve(Problem(Minimize(Bool())), solver = "ECOS"))
  expect_error(solve(Problem(Minimize(lambda_max(a))), solver = "ECOS"))
  expect_error(solve(Problem(Minimize(a)), solver = "SCS"))
})

test_that("Tests problems with reshape_expr", {
  skip_on_cran()
  # Test on scalars
  expect_equal(value(reshape_expr(1,c(1,1))), matrix(1))

  # Test vector to matrix
  x <- Variable(4)
  mat <- cbind(c(1,-1), c(2,-2))
  vec <- matrix(1:4)
  vec_mat <- cbind(c(1,2), c(3,4))
  expr <- reshape_expr(x,c(2,2))
  obj <- Minimize(sum_entries(mat %*% expr))
  prob <- Problem(obj, list(x == vec))
  result <- solve(prob)
  expect_equal(result$value, sum(mat %*% vec_mat))

  # Test on matrix to vector
  c <- 1:4
  expr <- reshape_expr(A,c(4,1))
  obj <- Minimize(t(expr) %*% c)
  constraints <- list(A == cbind(c(-1,-2), c(3,4)))
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  expect_equal(result$value, 20, tolerance = TOL)
  expect_equal(result$getValue(expr), matrix(c(-1,-2,3,4)))
  expect_equal(result$getValue(reshape_expr(expr,c(2,2))), cbind(c(-1,-2), c(3,4)))

  # Test matrix to matrix
  expr <- reshape_expr(C,c(2,3))
  mat <- rbind(c(1,-1), c(2,-2))
  C_mat <- rbind(c(1,4), c(2,5), c(3,6))
  obj <- Minimize(sum_entries(mat %*% expr))
  prob <- Problem(obj, list(C == C_mat))
  result <- solve(prob)
  reshaped = matrix(C_mat, nrow = 2, ncol = 3, byrow = FALSE)
  expect_equal(result$value, sum(mat %*% reshaped), tolerance = TOL)
  expect_equal(result$getValue(expr), reshaped, tolerance = TOL)

  # Test promoted expressions
  c <- cbind(c(1,-1), c(2,-2))
  expr <- reshape_expr(c * a,c(1,4))
  obj <- Minimize(expr %*% (1:4))
  prob <- Problem(obj, list(a == 2))
  result <- solve(prob)
  expect_equal(result$value, -6, tolerance = TOL)
  expect_equal(result$getValue(expr), 2*matrix(c, nrow = 1), tolerance = TOL)

  expr <- reshape_expr(c * a,c(4,1))
  obj <- Minimize(t(expr) %*% (1:4))
  prob <- Problem(obj, list(a == 2))
  result <- solve(prob)
  expect_equal(result$value, -6, tolerance = TOL)
  expect_equal(result$getValue(expr), 2*matrix(c, ncol = 1), tolerance = TOL)
})

test_that("Tests problems with vec", {
  skip_on_cran()
  c <- 1:4
  expr <- vec(A)
  obj <- Minimize(t(expr) %*% c)
  constraints <- list(A == cbind(c(-1,-2), c(3,4)))
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  expect_equal(result$value, 20, tolerance = TOL)
  expect_equal(result$getValue(expr), matrix(c(-1,-2,3,4)))
})

test_that("Test a problem with diag", {
  skip_on_cran()
  C <- Variable(3,3)
  obj <- Maximize(C[1,3])
  constraints <- list(diag(C) == 1, C[1,2] == 0.6, C[2,3] == -0.3, C == Variable(3, 3, PSD = TRUE))
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  expect_equal(result$value, 0.583151, tolerance = 1e-2)
})

# test_that("Test presolve parameters", {
#   # Test with parameters
#   gamma <- Parameter(nonneg = TRUE)
#   x <- Variable()
#   obj <- Minimize(x)
#   value(gamma) <- 0
#   prob <- Problem(obj, list(gamma == 1, x >= 0))
#   result <- solve(prob, solver = "SCS")
#   expect_equal(result$status, "infeasible")
#
#   value(gamma) <- 1
#   prob <- Problem(obj, list(gamma == 1, x >= 0))
#   result <- solve(prob, solver = "SCS")
#   expect_equal(result$status, "optimal")
# })

# test_that("Test that expressions with parameters are updated properly", {
#   x <- Variable()
#   y <- Variable()
#   x0 <- Parameter()
#
#   # Initial guess for x
#   value(x0) <- 2
#
#   # Make the constraints x^2 - y == 0
#   xSquared <- x0*x0 + 2*x0*(x-x0)
#   g <- xSquared - y
#
#   # Set up the problem
#   obj <- abs(x - 1)
#   prob <- Problem(Minimize(obj), list(g == 0))
#   result <- solve(prob)
#
#   value(x0) <- 1
#   xSquared <- x0*x0 + 2*x0*(x-x0)
#   g <- xSquared - y
#   prob <- Problem(Minimize(obj), list(g == 0))
#   result <- solve(prob)
#   # expect_equal(result$getValue(g), 0, tolerance = TOL)
#
#   # Test multiplication
#   prob <- Problem(Minimize(x0*x), list(x == 1))
#   value(x0) <- 2
#   result <- solve(prob)
#   x0@value <- 1
#   result <- solve(prob)
#   expect_equal(result$value, 1, tolerance = TOL)
# })

test_that("Test interaction of caching with changing constraints", {
  skip_on_cran()
  prob <- Problem(Minimize(a), list(a == 2, a >= 1))
  result <- solve(prob)
  expect_equal(result$value, 2, tolerance = TOL)

  prob@constraints[[1]] = (a == 1)
  result <- solve(prob)
  expect_equal(result$value, 1, tolerance = TOL)
})

test_that("Test positive definite constraints", {
  skip_on_cran()
  C <- Variable(3,3)
  obj <- Maximize(C[1,3])
  constraints <- list(diag(C) == 1, C[1,2] == 0.6, C[2,3] == -0.3, C == t(C), C %>>% 0)
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  expect_equal(result$value, 0.583151, tolerance = 1e-2)

  C <- Variable(2,2)
  obj <- Maximize(C[1,2])
  constraints <- list(C == 1, C %>>% rbind(c(2,0), c(0,2)))
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  expect_equal(result$status, "infeasible")

  C <- Variable(2, 2, symmetric = TRUE)
  obj <- Minimize(C[1,1])
  constraints <- list(C %<<% cbind(c(2,0), c(0,2)))
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  expect_equal(result$status, "unbounded")
})

test_that("Test the duals of PSD constraints", {
  skip_on_cran()
  # Test dual values with SCS
  C <- Variable(2, 2, symmetric = TRUE, name = "C")
  obj <- Maximize(C[1,1])
  constraints <- list(C %<<% cbind(c(2,0), c(0,2)))
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, 2, tolerance = 1e-4)

  psd_constr_dual <- result$getDualValue(constraints[[1]])
  C <- Variable(2, 2, symmetric = TRUE, name = "C")
  X <- Variable(2, 2, PSD = TRUE)
  obj <- Maximize(C[1,1])
  constraints <- list(X == cbind(c(2,0), c(0,2)) - C)
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$getDualValue(constraints[[1]]), psd_constr_dual, tolerance = 1e-3)

  # Test dual values with SCS that have off-diagonal entries
  C <- Variable(2, 2, symmetric = TRUE)
  obj <- Maximize(C[1,2] + C[2,1])
  constraints <- list(C %<<% cbind(c(2,0), c(0,2)), C >= 0)
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, 4, tolerance = 1e-3)

  psd_constr_dual <- result$getDualValue(constraints[[1]])
  C <- Variable(2, 2, symmetric = TRUE)
  X <- Variable(2, 2, PSD = TRUE)
  obj <- Maximize(C[1,2] + C[2,1])
  constraints <- list(X == cbind(c(2,0), c(0,2)) - C, C >= 0)
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$getDualValue(constraints[[1]]), psd_constr_dual, tolerance = 1e-3)
})

test_that("Test geo_mean function", {
  skip_on_cran()
    x <- Variable(2)
    cost <- geo_mean(x)
    prob <- Problem(Maximize(cost), list(x <= 1))
    result <- solve(prob)
    expect_equal(result$value, 1, tolerance = TOL)

    prob <- Problem(Maximize(cost), list(sum(x) <= 1))
    result <- solve(prob)
    expect_equal(result$getValue(x), matrix(c(0.5,0.5)), tolerance = TOL)

    x <- Variable(3,3)
    expect_error(geo_mean(x))

    x <- Variable(3,1)
    g <- geo_mean(x)
    ##expect_equal(g@w, gmp::as.bigq(1,3)*3, tolerance = TOL)

    x <- Variable(1,5)
    g <- geo_mean(x)
    ##expect_equal(g@w, gmp::as.bigq(1,5)*5, tolerance = TOL)

    ## Check that we get the right answer for max geo_mean(x) s.t. sum(x) <= 1
    p <- c(0.07, 0.12, 0.23, 0.19, 0.39)

    short_geo_mean <- function(x, p) {
        p <- as.numeric(p)/sum(p)
        x <- as.numeric(x)
        prod(x^p)
    }

    x <- Variable(5)
    prob <- Problem(Maximize(geo_mean(x, p)), list(sum(x) <= 1))
    result <- solve(prob)
    x <- as.numeric(result$getValue(x))
    x_true <- p/sum(p)

    expect_equal(result$value, result$getValue(geo_mean(x, p)), tolerance = TOL)
    expect_equal(result$value, short_geo_mean(x, p), tolerance = TOL)
    expect_equal(x, x_true, tolerance = 1e-3)

    ## Check that we get the right answer for max geo_mean(x) s.t. p_norm(x) <= 1
    x <- Variable(5)
    prob <- Problem(Maximize(geo_mean(x, p)), list(p_norm(x) <= 1))
    result <- solve(prob)
    x <- as.numeric(result$getValue(x))
    x_true <- sqrt(p/sum(p))

    expect_equal(result$value, result$getValue(geo_mean(x, p)), tolerance = 1e-3)
    expect_equal(result$value, short_geo_mean(x, p), tolerance = 1e-3)
    expect_equal(x, x_true, tolerance = 1e-3)

    ## The following 3 tests check vstack and hstack input to geo_mean
    ## The following 3 formulations should be equivalent
    n <- 5
    x_true <- rep(1,n)
    x <- Variable(n)

    result <- solve(Problem(Maximize(geo_mean(x)), list(x <= 1)))
    xval <- as.numeric(result$getValue(x))
    expect_equal(xval, x_true, tolerance = 1e-3)

    ##args <- list()
    ##for(i in 1:n)
    ##    args <- c(args, x[i])
    args <- lapply(seq_len(n), function(ind) x[ind])
    y <- do.call(vstack, args)
    result <- solve(Problem(Maximize(geo_mean(y)), list(x <= 1)))
    xval <- as.numeric(result$getValue(x))
    expect_equal(xval, x_true, tolerance = 1e-3)

    y <- do.call(hstack, args)
    result <- solve(Problem(Maximize(geo_mean(y)), list(x <= 1)))
    xval <- as.numeric(result$getValue(x))
    expect_equal(xval, x_true, tolerance = 1e-3)
})

test_that("Test p_norm function", {
  skip_on_cran()
    x <- Variable(3, name = "x")
    avec <- c(1.0, 2, 3)

    ## TODO: Add -1, 0.5, 0.3, -2.3 and testing positivity constraints

  for(p in c(1, 1.6, 1.3, 2, 1.99, 3, 3.7, Inf)) {
    prob <- Problem(Minimize(p_norm(x, p = p)), list(t(x) %*% avec >= 1))
    result <- solve(prob)

    # Formula is true for any a >= 0 with p > 1
    if(p == Inf)
      x_true <- rep(1, length(avec))/sum(avec)
    else if(p == 1) {
      # Only works for the particular a = c(1,2,3)
      x_true <- c(0,0,1.0/3)
    } else
      x_true <- avec^(1.0/(p-1)) / as.numeric(avec %*% (avec^(1.0/(p-1))))

    x_alg <- as.vector(result$getValue(x))
    expect_equal(x_alg, x_true, tolerance = 1e-3)
    expect_equal(result$value, .p_norm(x_alg, p), tolerance = TOL)
    expect_equal(.p_norm(x_alg, p), result$getValue(p_norm(x_alg, p)))
  }
})

test_that("Test p_norm concave", {
  skip_on_cran()
  x <- Variable(3, name = "x")

  # Test positivity constraints
  a <- c(-1.0, 2, 3)
  for(p in c(-1, 0.5, 0.3, -2.3)) {
    prob <- Problem(Minimize(sum_entries(abs(x-a))), list(p_norm(x,p) >= 0))
    result <- solve(prob)
    expect_equal(result$value, 1, tolerance = TOL)
  }

  a <- c(1.0, 2, 3)
  for(p in c(-1, 0.5, 0.3, -2.3)) {
    prob <- Problem(Minimize(sum_entries(abs(x-a))), list(p_norm(x,p) >= 0))
    result <- solve(prob)
    expect_equal(result$value, 0, tolerance = TOL)
  }
})

test_that("Test power function", {
  skip_on_cran()
  x <- Variable()
  prob <- Problem(Minimize(power(x, 1.7) + power(x, -2.3) - power(x, 0.45)))
  result <- solve(prob)
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
  result <- solve(prob)
  expect_equal(result$value, 8, tolerance = TOL)
})
