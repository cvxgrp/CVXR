## DPP (Disciplined Parameterized Programming) parity tests
## Ported from CVXPY: cvxpy/tests/test_dpp.py (TestDcp class)
## Tests cover is_dpp detection, chain inspection, and paper examples.

# ══════════════════════════════════════════════════════════════════
# CRITICAL: is_dpp detection tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_dpp.py::TestDcp::test_multiply_param_and_variable_is_dpp
test_that("DPP parity: multiply param and variable is DPP", {
  ## CVXPY: test_multiply_param_and_variable_is_dpp
  x <- Parameter()
  y <- Variable()
  product <- x * y
  expect_true(is_dpp(product))
  expect_true(is_dcp(product))
})

## @cvxpy test_dpp.py::TestDcp::test_multiply_variable_and_param_is_dpp
test_that("DPP parity: multiply variable and param is DPP", {
  ## CVXPY: test_multiply_variable_and_param_is_dpp
  ## cp.multiply(y, x) -> y * x (elementwise)
  x <- Parameter()
  y <- Variable()
  product <- y * x
  expect_true(is_dpp(product))
  expect_true(is_dcp(product))
})

## @cvxpy test_dpp.py::TestDcp::test_multiply_nonlinear_param_and_variable_is_not_dpp
test_that("DPP parity: multiply nonlinear param and variable is NOT DPP", {
  ## CVXPY: test_multiply_nonlinear_param_and_variable_is_not_dpp
  x <- Parameter()
  y <- Variable()
  product <- exp(x) * y
  expect_false(is_dpp(product))
})

## @cvxpy test_dpp.py::TestDcp::test_multiply_affine_param_and_variable_is_dpp
test_that("DPP parity: multiply affine param and variable is DPP", {
  ## CVXPY: test_multiply_affine_param_and_variable_is_dpp
  x <- Parameter()
  y <- Variable()
  product <- (x + x) * y
  expect_true(is_dpp(product))
  expect_true(is_dcp(product))
})

## @cvxpy test_dpp.py::TestDcp::test_solve_dpp_problem
test_that("DPP parity: solve DPP problem and re-solve with new param", {
  ## CVXPY: test_solve_dpp_problem
  x <- Parameter()
  value(x) <- 5
  y <- Variable()
  problem <- Problem(Minimize(x + y), list(x == y))
  expect_true(is_dpp(problem))
  expect_true(is_dcp(problem))
  result1 <- psolve(problem, solver = "SCS")
  expect_equal(result1, 10, tolerance = 1e-3)

  value(x) <- 3
  result2 <- psolve(problem, solver = "SCS")
  expect_equal(result2, 6, tolerance = 1e-3)
})

## @cvxpy test_dpp.py::TestDcp::test_chain_data_for_dpp_problem_does_not_eval_params
test_that("DPP parity: chain for DPP problem does not contain EvalParams", {
  ## CVXPY: test_chain_data_for_dpp_problem_does_not_eval_params
  x <- Parameter()
  value(x) <- 5
  y <- Variable()
  problem <- Problem(Minimize(x + y), list(x == y))

  pd <- problem_data(problem, solver = "SCS")
  chain <- pd$chain
  has_eval_params <- any(vapply(chain@reductions, function(r) {
    S7_inherits(r, EvalParams)
  }, logical(1L)))
  expect_false(has_eval_params)
})

## @cvxpy test_dpp.py::TestDcp::test_paper_example_logreg_is_dpp
test_that("DPP parity: paper example (logreg) is DPP", {
  ## CVXPY: test_paper_example_logreg_is_dpp
  N <- 3L; n <- 2L
  beta <- Variable(c(n, 1L))
  b <- Variable(c(1L, 1L))
  X <- Parameter(c(N, n))
  Y <- matrix(1, N, 1)
  lambd1 <- Parameter(nonneg = TRUE)
  lambd2 <- Parameter(nonneg = TRUE)

  ## log_likelihood = (1/N) * sum(Y * (X @ beta + b) -
  ##   log_sum_exp(hstack(zeros, X @ beta + b).T, axis=0, keepdims=TRUE).T)
  Xbeta_b <- X %*% beta + b
  zeros_N1 <- matrix(0, N, 1)
  stacked <- hstack(zeros_N1, Xbeta_b)  # (N, 2)
  ## CVXPY: .T then axis=0, keepdims=TRUE then .T
  ## In R: transpose -> (2, N), log_sum_exp(axis=0, keepdims=TRUE) -> (1, N), transpose -> (N, 1)
  lse_part <- t(log_sum_exp(t(stacked), axis = 1L, keepdims = TRUE))
  log_likelihood <- (1.0 / N) * sum_entries(Y * Xbeta_b - lse_part)
  regularization <- -lambd1 * norm1(beta) - lambd2 * sum_squares(beta)
  problem <- Problem(Maximize(log_likelihood + regularization))

  expect_true(is_dpp(log_likelihood))
  expect_true(is_dcp(problem))
  expect_true(is_dpp(problem))
})

## @cvxpy test_dpp.py::TestDcp::test_paper_example_stoch_control
test_that("DPP parity: paper example (stochastic control) is DPP", {
  ## CVXPY: test_paper_example_stoch_control
  n <- 3L; m <- 3L
  x <- Parameter(c(n, 1L))
  P_sqrt <- Parameter(c(m, m))
  P_21 <- Parameter(c(n, m))
  q <- Parameter(c(m, 1L))
  u <- Variable(c(m, 1L))
  y <- Variable(c(n, 1L))
  objective <- 0.5 * sum_squares(P_sqrt %*% u) + t(x) %*% y + t(q) %*% u
  problem <- Problem(Minimize(objective),
                     list(cvxr_norm(u) <= 0.5, y == P_21 %*% u))
  expect_true(is_dpp(problem))
  expect_true(is_dcp(problem))
})

## @cvxpy test_dpp.py::TestDcp::test_paper_example_relu
test_that("DPP parity: paper example (ReLU) solves correctly", {
  ## CVXPY: test_paper_example_relu
  n <- 2L
  x <- Parameter(n)
  y <- Variable(n)
  objective <- Minimize(sum_squares(y - x))
  constraints <- list(y >= 0)
  problem <- Problem(objective, constraints)
  expect_true(is_dpp(problem))

  ## When x = [5, 5], ReLU(x) = [5, 5]
  value(x) <- c(5, 5)
  psolve(problem, solver = "SCS")
  expect_equal(as.numeric(value(y)), c(5, 5), tolerance = 1e-4)

  ## When x = [-4, -4], ReLU(x) = [0, 0]
  value(x) <- c(-4, -4)
  psolve(problem, solver = "SCS")
  expect_equal(as.numeric(value(y)), c(0, 0), tolerance = 1e-4)
})

## @cvxpy test_dpp.py::TestDcp::test_paper_example_opt_net_qp
test_that("DPP parity: paper example (OptNet QP) is DPP", {
  ## CVXPY: test_paper_example_opt_net_qp
  m <- 3L; n <- 2L
  G <- Parameter(c(m, n))
  h <- Parameter(c(m, 1L))
  p <- Parameter(c(n, 1L))
  y <- Variable(c(n, 1L))
  objective <- Minimize(0.5 * sum_squares(y - p))
  constraints <- list(G %*% y <= h)
  problem <- Problem(objective, constraints)
  expect_true(is_dpp(problem))
})

## @cvxpy test_dpp.py::TestDcp::test_paper_example_ellipsoidal_constraints
test_that("DPP parity: paper example (ellipsoidal constraints) is DPP", {
  ## CVXPY: test_paper_example_ellipsoidal_constraints
  n <- 2L
  A_sqrt <- Parameter(c(n, n))
  z <- Parameter(n)
  p <- Parameter(n)
  y <- Variable(n)
  slack <- Variable(c(n, 1L))
  objective <- Minimize(0.5 * sum_squares(y - p))
  constraints <- list(
    0.5 * sum_squares(A_sqrt %*% slack) <= 1,
    slack == y - z
  )
  problem <- Problem(objective, constraints)
  expect_true(is_dpp(problem))
})


# ══════════════════════════════════════════════════════════════════
# MEDIUM: is_dpp detection tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_dpp.py::TestDcp::test_multiply_scalar_params_not_dpp
test_that("DPP parity: multiply scalar params is NOT DPP", {
  ## CVXPY: test_multiply_scalar_params_not_dpp
  x <- Parameter()
  product <- x * x
  expect_false(is_dpp(product))
  expect_true(is_dcp(product))
})

## @cvxpy test_dpp.py::TestDcp::test_matmul_params_not_dpp
test_that("DPP parity: matmul params is NOT DPP", {
  ## CVXPY: test_matmul_params_not_dpp
  X <- Parameter(c(4L, 4L))
  product <- X %*% X
  expect_true(is_dcp(product))
  expect_false(is_dpp(product))
})

## @cvxpy test_dpp.py::TestDcp::test_nonlinear_equality_not_dpp
test_that("DPP parity: nonlinear equality constraint is NOT DPP", {
  ## CVXPY: test_nonlinear_equality_not_dpp
  x <- Variable()
  a <- Parameter()
  constraint <- list(x == cvxr_norm(a))
  ## In CVXPY: constraint[0].is_dcp(dpp=True) -> False
  ## In CVXR: is_dpp(constraint) uses with_dpp_scope(is_dcp(constraint))
  expect_false(is_dpp(constraint[[1]]))
  problem <- Problem(Minimize(0), constraint)
  expect_false(is_dpp(problem))
})

## @cvxpy test_dpp.py::TestDcp::test_nonconvex_inequality_not_dpp
test_that("DPP parity: nonconvex inequality constraint is NOT DPP", {
  ## CVXPY: test_nonconvex_inequality_not_dpp
  x <- Variable()
  a <- Parameter()
  constraint <- list(x <= cvxr_norm(a))
  ## x <= norm(a): in DPP scope, norm(a) is convex in a (param treated as variable),
  ## so the constraint "x - norm(a) <= 0" has nonconvex part on RHS
  expect_false(is_dpp(constraint[[1]]))
  problem <- Problem(Minimize(0), constraint)
  expect_false(is_dpp(problem))
})

## @cvxpy test_dpp.py::TestDcp::test_non_dcp_expression_is_not_dpp
test_that("DPP parity: non-DCP expression is NOT DPP", {
  ## CVXPY: test_non_dcp_expression_is_not_dpp
  x <- Parameter()
  expr <- exp(log(x))
  expect_false(is_dpp(expr))
})

## @cvxpy test_dpp.py::TestDcp::test_can_solve_non_dpp_problem
test_that("DPP parity: can solve non-DPP problem", {
  ## CVXPY: test_can_solve_non_dpp_problem
  x <- Parameter()
  value(x) <- 5
  y <- Variable()
  problem <- Problem(Minimize(x * x), list(x == y))
  expect_false(is_dpp(problem))
  expect_true(is_dcp(problem))

  ## Should still solve (via EvalParams path)
  result1 <- suppressWarnings(psolve(problem, solver = "SCS"))
  expect_equal(result1, 25, tolerance = 1e-3)

  value(x) <- 3
  result2 <- suppressWarnings(psolve(problem, solver = "SCS"))
  expect_equal(result2, 9, tolerance = 1e-3)
})

## @cvxpy test_dpp.py::TestDcp::test_param_quad_form_not_dpp
test_that("DPP parity: param quad_form is NOT DPP", {
  ## CVXPY: test_param_quad_form_not_dpp
  x <- Variable(c(2L, 1L))
  P <- Parameter(c(2L, 2L), PSD = TRUE)
  value(P) <- diag(2)
  qf <- quad_form(x, P)
  expect_false(is_dpp(qf))
  expect_true(is_dcp(qf))
})

## @cvxpy test_dpp.py::TestDcp::test_const_quad_form_is_dpp
test_that("DPP parity: constant quad_form IS DPP", {
  ## CVXPY: test_const_quad_form_is_dpp
  x <- Variable(c(2L, 1L))
  P <- diag(2)
  qf <- quad_form(x, P)
  expect_true(is_dpp(qf))
  expect_true(is_dcp(qf))
})

## @cvxpy test_dpp.py::TestDcp::test_non_dpp_powers
test_that("DPP parity: non-DPP powers solve correctly", {
  ## CVXPY: test_non_dpp_powers
  ## Case 1: s in objective (affine in param, so DPP)
  s <- Parameter(1L, nonneg = TRUE)
  x <- Variable(1L)
  obj <- Maximize(x + s)
  cons <- list(x <= 1)
  prob <- Problem(obj, cons)
  value(s) <- 1.0
  suppressWarnings(psolve(prob, solver = "SCS"))
  expect_equal(value(prob), 2.0, tolerance = 1e-3)

  ## Case 2: s^2 in objective (nonlinear in param, NOT DPP)
  s <- Parameter(1L, nonneg = TRUE)
  x <- Variable(1L)
  obj <- Maximize(x + power(s, 2))
  cons <- list(x <= 1)
  prob <- Problem(obj, cons)
  value(s) <- 1.0
  suppressWarnings(psolve(prob, solver = "SCS"))
  expect_equal(value(prob), 2.0, tolerance = 1e-3)

  ## Case 3: multiply(x, s^2) (nonlinear param * variable, NOT DPP)
  s <- Parameter(1L, nonneg = TRUE)
  x <- Variable(1L)
  obj <- Maximize(x * power(s, 2))
  cons <- list(x <= 1)
  prob <- Problem(obj, cons)
  value(s) <- 1.0
  suppressWarnings(psolve(prob, solver = "SCS"))
  expect_equal(value(prob), 1.0, tolerance = 1e-3)
})

## @cvxpy test_dpp.py::TestDgp::test_param_monomial_is_dpp
test_that("DPP parity: param monomial is DPP (DGP context)", {
  ## CVXPY: test_param_monomial_is_dpp
  alpha <- Parameter(pos = TRUE)
  beta <- Parameter(pos = TRUE)
  kappa <- Parameter(pos = TRUE)

  monomial <- power(alpha, 1.2) * power(beta, 0.5) * power(kappa, 3) * power(kappa, 2)
  ## In DGP DPP context: this should be DGP-DPP compliant
  ## is_dgp(dpp=True) in CVXPY -> .is_dgp_dpp in CVXR (checks log-log convexity in DPP scope)
  ## We check that the expression is log-log convex in DPP scope
  expect_true(with_dpp_scope(is_log_log_convex(monomial)))
})

## @cvxpy test_dpp.py::TestDgp::test_param_posynomial_is_dpp
test_that("DPP parity: param posynomial is DPP (DGP context)", {
  ## CVXPY: test_param_posynomial_is_dpp
  alpha <- Parameter(pos = TRUE)
  beta <- Parameter(pos = TRUE)
  kappa <- Parameter(pos = TRUE)

  monomial <- power(alpha, 1.2) * power(beta, 0.5) * power(kappa, 3) * power(kappa, 2)
  posynomial <- monomial + power(alpha, 2) * power(beta, 3)
  ## In DGP DPP context: posynomial of positive params is DGP-DPP
  expect_true(with_dpp_scope(is_log_log_convex(posynomial)))
})

# ══════════════════════════════════════════════════════════════════
# ADDITIONAL: chain inspection and chain data tests
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_dpp.py::TestDcp::test_chain_data_for_non_dpp_problem_evals_params
test_that("DPP parity: chain for non-DPP problem DOES contain EvalParams", {
  ## CVXPY: test_chain_data_for_non_dpp_problem_evals_params
  x <- Parameter()
  value(x) <- 5
  y <- Variable()
  problem <- Problem(Minimize(x * x), list(x == y))

  pd <- suppressWarnings(problem_data(problem, solver = "SCS"))
  chain <- pd$chain
  expect_false(is_dpp(problem))

  has_eval_params <- any(vapply(chain@reductions, function(r) {
    S7_inherits(r, EvalParams)
  }, logical(1L)))
  expect_true(has_eval_params)
})

## @cvxpy test_dpp.py::TestDcp::test_paper_example_is_dpp
test_that("DPP parity: paper example (Fx-g + lambda*x) is DPP", {
  ## CVXPY: test_paper_example_is_dpp
  F_param <- Parameter(c(2L, 2L))
  x <- Variable(c(2L, 1L))
  g <- Parameter(c(2L, 1L))
  lambd <- Parameter(nonneg = TRUE)
  objective <- cvxr_norm(F_param %*% x - g) + lambd * cvxr_norm(x)
  constraints <- list(x >= 0)
  problem <- Problem(Minimize(objective), constraints)
  expect_true(is_dpp(objective))
  expect_true(is_dpp(constraints[[1]]))
  expect_true(is_dpp(problem))
})

## @cvxpy test_dpp.py::TestDcp::test_multiply_param_plus_var_times_const
test_that("DPP parity: multiply param+var times const is DPP", {
  ## CVXPY: test_multiply_param_plus_var_times_const
  x <- Parameter()
  y <- Variable()
  product <- (x + y) * 5
  expect_true(is_convex(product))
  expect_true(is_dcp(product))
  expect_true(is_dpp(product))
})

## @cvxpy test_dpp.py::TestDcp::test_solve_multiply_param_plus_var_times_const
test_that("DPP parity: solve (param+var)*const problem", {
  ## CVXPY: test_solve_multiply_param_plus_var_times_const
  x <- Parameter()
  y <- Variable()
  product <- (x + y) * 5
  expect_true(is_dpp(product))
  value(x) <- 2.0
  problem <- Problem(Minimize(product), list(y == 1))
  result <- psolve(problem, solver = "SCS")
  expect_equal(result, 15, tolerance = 1e-3)
})

## @cvxpy test_dpp.py::TestDcp::test_multiply_nonlinear_nonneg_param_and_nonneg_variable_is_not_dpp
test_that("DPP parity: multiply nonlinear nonneg param and nonneg variable is NOT DPP", {
  ## CVXPY: test_multiply_nonlinear_nonneg_param_and_nonneg_variable_is_not_dpp
  x <- Parameter(nonneg = TRUE)
  y <- Variable(nonneg = TRUE)
  product <- exp(x) * y
  expect_false(is_dpp(product))
  expect_true(is_dcp(product))
})

## @cvxpy test_dpp.py::TestDcp::test_multiply_param_and_nonlinear_variable_is_dpp
test_that("DPP parity: multiply param and nonlinear variable IS DPP", {
  ## CVXPY: test_multiply_param_and_nonlinear_variable_is_dpp
  x <- Parameter(nonneg = TRUE)
  y <- Variable()
  product <- x * exp(y)
  expect_true(is_convex(product))
  expect_true(is_dcp(product))
  expect_true(is_dpp(product))
})

# ══════════════════════════════════════════════════════════════════
# DPP solve: paper examples with actual numerical verification
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("DPP parity: stochastic control solves correctly with parameter changes", {
  ## Extended from test_paper_example_stoch_control: actually solve the problem
  n <- 3L; m <- 3L
  x_param <- Parameter(c(n, 1L))
  P_sqrt <- Parameter(c(m, m))
  P_21 <- Parameter(c(n, m))
  q_param <- Parameter(c(m, 1L))
  u <- Variable(c(m, 1L))
  y_var <- Variable(c(n, 1L))

  objective <- 0.5 * sum_squares(P_sqrt %*% u) + t(x_param) %*% y_var + t(q_param) %*% u
  problem <- Problem(Minimize(objective),
                     list(cvxr_norm(u) <= 0.5, y_var == P_21 %*% u))
  expect_true(is_dpp(problem))

  ## Set parameter values and solve
  value(x_param) <- matrix(c(1, 0, 0), 3, 1)
  value(P_sqrt) <- diag(3)
  value(P_21) <- diag(3)
  value(q_param) <- matrix(0, 3, 1)

  result <- psolve(problem, solver = "SCS")
  expect_true(is.finite(result))

  ## Re-solve with different parameters (DPP fast path)
  value(x_param) <- matrix(c(0, 1, 0), 3, 1)
  result2 <- psolve(problem, solver = "SCS")
  expect_true(is.finite(result2))
})

## @cvxpy NONE
test_that("DPP parity: OptNet QP solves correctly with parameter changes", {
  ## Extended from test_paper_example_opt_net_qp: actually solve
  m <- 3L; n <- 2L
  G <- Parameter(c(m, n))
  h <- Parameter(c(m, 1L))
  p_param <- Parameter(c(n, 1L))
  y <- Variable(c(n, 1L))
  objective <- Minimize(0.5 * sum_squares(y - p_param))
  constraints <- list(G %*% y <= h)
  problem <- Problem(objective, constraints)
  expect_true(is_dpp(problem))

  ## Set params: G = identity-ish, h = [1,1,1], p = [0.5, 0.5]
  value(G) <- rbind(diag(2), c(0, 0))
  value(h) <- matrix(c(1, 1, 1), 3, 1)
  value(p_param) <- matrix(c(0.5, 0.5), 2, 1)
  result1 <- psolve(problem, solver = "SCS")
  ## Unconstrained min is at y = p = [0.5, 0.5], which satisfies G*y <= h
  expect_equal(result1, 0.0, tolerance = 1e-3)

  ## Change p to be outside the feasible region
  value(p_param) <- matrix(c(2, 2), 2, 1)
  result2 <- psolve(problem, solver = "SCS")
  ## y is constrained by Gy <= h, so optimal is clamped
  expect_true(result2 > 0)  # can't reach p=[2,2]
})

## @cvxpy NONE
test_that("DPP parity: ellipsoidal constraints solve with parameter changes", {
  ## Extended from test_paper_example_ellipsoidal_constraints
  n <- 2L
  A_sqrt <- Parameter(c(n, n))
  z <- Parameter(n)
  p_param <- Parameter(n)
  y <- Variable(n)
  slack <- Variable(c(n, 1L))
  objective <- Minimize(0.5 * sum_squares(y - p_param))
  constraints <- list(
    0.5 * sum_squares(A_sqrt %*% slack) <= 1,
    slack == y - z
  )
  problem <- Problem(objective, constraints)
  expect_true(is_dpp(problem))

  ## Solve with identity A_sqrt, z=0, p=[3,3]
  value(A_sqrt) <- diag(2)
  value(z) <- c(0, 0)
  value(p_param) <- c(3, 3)
  result1 <- psolve(problem, solver = "CLARABEL")
  expect_true(is.finite(result1))

  ## Re-solve with different parameters
  value(p_param) <- c(0, 0)
  result2 <- psolve(problem, solver = "CLARABEL")
  expect_equal(result2, 0.0, tolerance = 1e-3)
})

# ══════════════════════════════════════════════════════════════════
# DPP fast path: verify caching behavior
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("DPP parity: DPP problem caches param_prog after first solve", {
  x <- Parameter()
  value(x) <- 5
  y <- Variable()
  problem <- Problem(Minimize(x + y), list(x == y))
  expect_true(is_dpp(problem))

  ## First solve caches
  psolve(problem, solver = "SCS")
  expect_false(is.null(problem@.cache$param_prog))

  ## Second solve with different params uses fast path
  value(x) <- 10
  result <- psolve(problem, solver = "SCS")
  expect_equal(result, 20, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("DPP parity: non-DPP problem falls back to EvalParams (no caching)", {
  x <- Parameter()
  value(x) <- 2
  y <- Variable()
  ## x * x is not DPP (param * param)
  problem <- Problem(Minimize(x * x), list(x == y))
  expect_false(is_dpp(problem))

  suppressWarnings(psolve(problem, solver = "SCS"))
  ## Non-DPP: no param_prog cached
  expect_null(problem@.cache$param_prog)
})

# ══════════════════════════════════════════════════════════════════
# DPP with QP problems
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("DPP parity: QP with parametric linear term is DPP", {
  ## Lasso-like problem: min ||Ax - b||^2 + gamma * ||x||_1
  set.seed(42)
  n <- 5L; m <- 3L
  A <- matrix(rnorm(m * n), m, n)
  b <- rnorm(m)
  x <- Variable(n)
  gamma <- Parameter(nonneg = TRUE)

  prob <- Problem(Minimize(sum_squares(A %*% x - b) + gamma * p_norm(x, 1)))
  expect_true(is_dpp(prob))

  ## Solve with gamma = 0 (least squares)
  value(gamma) <- 0
  val0 <- psolve(prob, solver = "CLARABEL")
  expect_true(is.finite(val0))

  ## Solve with gamma = 1 (lasso)
  value(gamma) <- 1
  val1 <- psolve(prob, solver = "CLARABEL")
  expect_true(val1 >= val0)  # more regularization -> larger objective

  ## Solve with gamma = 10
  value(gamma) <- 10
  val10 <- psolve(prob, solver = "CLARABEL")
  expect_true(val10 >= val1)
})

## @cvxpy NONE
test_that("DPP parity: QP parametric objective and constraints", {
  ## Ridge regression with parametric bound
  n <- 3L
  x <- Variable(n)
  c_param <- Parameter(n)
  bound <- Parameter(nonneg = TRUE)

  prob <- Problem(Minimize(sum_squares(x) + t(c_param) %*% x),
                  list(sum_entries(x) >= bound))
  expect_true(is_dpp(prob))

  value(c_param) <- c(-2, -1, 0)
  value(bound) <- 1
  val1 <- psolve(prob, solver = "CLARABEL")
  expect_true(is.finite(val1))

  ## Change parameters (uses DPP fast path)
  value(c_param) <- c(0, 0, 0)
  value(bound) <- 0.5
  val2 <- psolve(prob, solver = "CLARABEL")
  expect_true(is.finite(val2))
})

# ══════════════════════════════════════════════════════════════════
# DPP with multiple solvers
# ══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("DPP parity: DPP re-solve produces consistent results across solvers", {
  ## Solve the same DPP problem with Clarabel and SCS
  p <- Parameter()
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= p))
  expect_true(is_dpp(prob))

  value(p) <- 3
  val_scs <- psolve(prob, solver = "SCS")
  expect_equal(val_scs, 3.0, tolerance = 1e-3)

  ## Now create a fresh problem (to avoid cache interference) and solve with Clarabel
  p2 <- Parameter()
  x2 <- Variable()
  prob2 <- Problem(Minimize(x2), list(x2 >= p2))

  value(p2) <- 3
  val_clar <- psolve(prob2, solver = "CLARABEL")
  expect_equal(val_clar, 3.0, tolerance = 1e-3)

  ## Re-solve with different param (DPP fast path for both)
  value(p) <- 7
  val_scs2 <- psolve(prob, solver = "SCS")
  expect_equal(val_scs2, 7.0, tolerance = 1e-3)

  value(p2) <- 7
  val_clar2 <- psolve(prob2, solver = "CLARABEL")
  expect_equal(val_clar2, 7.0, tolerance = 1e-3)
})
