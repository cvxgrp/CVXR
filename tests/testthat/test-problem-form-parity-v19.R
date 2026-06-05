## CVXPY 1.9 ProblemForm / solving-chain parity.
##
## CVXR does not expose CVXPY's ProblemForm class directly. These tests cover
## the equivalent solver-chain predicates, cone detection/stuffing, and GP
## dispatch behavior that CVXPY's test_problem_form.py exercises.

problem_form_cone_names <- function(prob) {
  vapply(CVXR:::.required_cone_types(prob), function(cls) cls@name, character(1L))
}

problem_form_chain_solver <- function(prob, solver = NULL, gp = FALSE) {
  chain <- CVXR:::construct_solving_chain(prob, solver = solver, gp = gp)
  CVXR:::solver_name(chain@solver)
}

problem_form_is_commercial_solver <- function(solver) {
  solver %in% c("MOSEK", "GUROBI", "CPLEX", "XPRESS")
}

problem_form_expect_solver <- function(actual, expected) {
  expect_true(actual == expected || problem_form_is_commercial_solver(actual),
              info = sprintf("expected %s or commercial solver, got %s",
                             expected, actual))
}

problem_form_is_quad_stuffing <- function(reduction) {
  S7::S7_inherits(reduction, ConeMatrixStuffing) && isTRUE(reduction@quad_obj)
}

## @cvxpy test_problem_form.py::TestProblemForm::test_lp
test_that("ProblemForm parity: LP structure has linear objective and base cones", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1, x[1] == 2))

  expect_false(CVXR:::has_quadratic_term(prob@objective@args[[1L]]))
  expect_false(is_mixed_integer(prob))
  expect_true(length(prob@constraints) > 0L)

  data <- problem_data(prob, solver = "CLARABEL")$data
  dims <- data[[CVXR:::SD_DIMS]]
  expect_equal(dims@zero, 1L)
  expect_equal(dims@nonneg, 2L)
  expect_length(problem_form_cone_names(prob), 0L)
})

## @cvxpy test_problem_form.py::TestProblemForm::test_qp
test_that("ProblemForm parity: QP structure detects quadratic objective", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x)), list(x >= 1))

  expect_true(CVXR:::has_quadratic_term(prob@objective@args[[1L]]))
  expect_true(length(prob@constraints) > 0L)
})

## @cvxpy test_problem_form.py::TestProblemForm::test_qp_filtering
test_that("ProblemForm parity: QP path filters quadratic SOCs but keeps other cones", {
  x <- Variable()
  qp <- Problem(Minimize(sum_squares(x)))
  qp_cone <- reduction_apply(CVXR:::Dcp2Cone(quad_obj = TRUE), qp)[[1L]]
  expect_false(any(vapply(qp_cone@constraints, S7::S7_inherits, logical(1L), SOC)))

  y <- Variable()
  exp_qp <- Problem(Minimize(quad_over_lin(exp(y), 1)))
  exp_cone <- reduction_apply(CVXR:::Dcp2Cone(quad_obj = TRUE), exp_qp)[[1L]]
  expect_true(any(vapply(exp_cone@constraints, S7::S7_inherits, logical(1L), ExpCone)))
  expect_false(any(vapply(exp_cone@constraints, S7::S7_inherits, logical(1L), SOC)))

  z <- Variable()
  nested_qp <- Problem(Minimize(sum_squares(sum_squares(z))))
  nested_cone <- reduction_apply(CVXR:::Dcp2Cone(quad_obj = TRUE), nested_qp)[[1L]]
  expect_true(any(vapply(nested_cone@constraints, S7::S7_inherits, logical(1L), SOC)))
})

## @cvxpy test_problem_form.py::TestProblemForm::test_quad_over_lin_constant_denom_qp_path
test_that("ProblemForm parity: quad_over_lin with constant denominator stays on QP path", {
  x <- Variable()
  prob <- Problem(Minimize(quad_over_lin(exp(x), 1)))
  conic_prob <- reduction_apply(CVXR:::Dcp2Cone(quad_obj = TRUE), prob)[[1L]]

  expect_true(any(vapply(conic_prob@constraints, S7::S7_inherits, logical(1L), ExpCone)))
  expect_false(any(vapply(conic_prob@constraints, S7::S7_inherits, logical(1L), SOC)))
})

## @cvxpy test_problem_form.py::TestProblemForm::test_socp
test_that("ProblemForm parity: SOCP structure reports SOC cones", {
  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 2)), list(x >= 1))

  expect_true("SOC" %in% problem_form_cone_names(prob))
  expect_true(length(prob@constraints) > 0L)
})

## @cvxpy test_problem_form.py::TestProblemForm::test_sdp
test_that("ProblemForm parity: SDP structure reports PSD cones", {
  X <- Variable(c(2, 2), symmetric = TRUE)
  prob <- Problem(Minimize(lambda_max(X)), list(X %>>% 0))

  expect_true("PSD" %in% problem_form_cone_names(prob))
  expect_true(length(prob@constraints) > 0L)
})

## @cvxpy test_problem_form.py::TestProblemForm::test_exp_cone
test_that("ProblemForm parity: exponential atoms report ExpCone", {
  x <- Variable(2, pos = TRUE)
  prob <- Problem(Maximize(sum_entries(log(x))), list(x <= 2))

  expect_true("ExpCone" %in% problem_form_cone_names(prob))
})

## @cvxpy test_problem_form.py::TestProblemForm::test_pow_cone
test_that("ProblemForm parity: explicit power-cone constraints report PowCone3D", {
  x <- Variable()
  y <- Variable()
  z <- Variable()
  prob <- Problem(Minimize(x), list(PowCone3D(x, y, z, 0.5)))

  expect_true("PowCone3D" %in% problem_form_cone_names(prob))
})

## @cvxpy test_problem_form.py::TestProblemForm::test_mixed_integer
test_that("ProblemForm parity: integer and boolean variables are mixed-integer", {
  x <- Variable(2, integer = TRUE)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 0, x <= 5))
  expect_true(is_mixed_integer(prob))

  y <- Variable(2, boolean = TRUE)
  prob2 <- Problem(Minimize(sum_entries(y)), list(y >= 0))
  expect_true(is_mixed_integer(prob2))
})

## @cvxpy test_problem_form.py::TestProblemForm::test_unconstrained_linear
test_that("ProblemForm parity: unconstrained linear objective has no cones", {
  x <- Variable()
  prob <- Problem(Minimize(x))

  expect_false(CVXR:::has_quadratic_term(prob@objective@args[[1L]]))
  expect_length(problem_form_cone_names(prob), 0L)
  expect_false(is_mixed_integer(prob))
  expect_length(prob@constraints, 0L)
})

## @cvxpy test_problem_form.py::TestProblemForm::test_sdp_nsd_variable
test_that("ProblemForm parity: NSD variables create PSD cone dimensions", {
  X <- Variable(c(2, 2), NSD = TRUE)
  prob <- Problem(Maximize(matrix_trace(X)))
  data <- problem_data(prob, solver = "SCS")$data

  expect_equal(data[[CVXR:::SD_DIMS]]@psd, 2L)
})

## @cvxpy test_problem_form.py::TestPickDefaultSolver::test_lp_gets_clarabel test_problem_form.py::TestPickDefaultSolver::test_qp_gets_osqp test_problem_form.py::TestPickDefaultSolver::test_sdp_gets_scs test_problem_form.py::TestPickDefaultSolver::test_socp_gets_clarabel test_problem_form.py::TestPickDefaultSolver::test_mi_gets_highs_or_premium
test_that("ProblemForm parity: default solver selection matches CVXPY policy modulo commercial solvers", {
  x <- Variable(2)
  lp <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  problem_form_expect_solver(problem_form_chain_solver(lp), "CLARABEL")

  q <- Variable(2)
  qp <- Problem(Minimize(sum_squares(q)), list(q >= 1))
  problem_form_expect_solver(problem_form_chain_solver(qp), "OSQP")

  X <- Variable(c(2, 2), symmetric = TRUE)
  sdp <- Problem(Minimize(lambda_max(X)), list(X %>>% 0))
  problem_form_expect_solver(problem_form_chain_solver(sdp), "SCS")

  s <- Variable(2)
  socp <- Problem(Minimize(p_norm(s, 2)), list(s >= 1))
  problem_form_expect_solver(problem_form_chain_solver(socp), "CLARABEL")

  m <- Variable(2, integer = TRUE)
  mip <- Problem(Minimize(sum_entries(m)), list(m >= 0, m <= 5))
  mi_solver <- problem_form_chain_solver(mip)
  expect_true(mi_solver == "HIGHS" || problem_form_is_commercial_solver(mi_solver),
              info = sprintf("expected HIGHS or commercial solver, got %s", mi_solver))
})

## @cvxpy test_problem_form.py::TestResolveAndBuildChain::test_solver_none_picks_default
test_that("ProblemForm parity: solver NULL builds a chain and solves LP", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  val <- psolve(prob, solver = NULL)
  expect_equal(val, 2, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-5)
})

## @cvxpy test_problem_form.py::TestResolveAndBuildChain::test_solver_string_clarabel
test_that("ProblemForm parity: named CLARABEL chain solves SOCP", {
  skip_if_not("CLARABEL" %in% installed_solvers())

  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 2)), list(x >= 1))
  expect_equal(problem_form_chain_solver(prob, solver = "CLARABEL"), "CLARABEL")

  val <- psolve(prob, solver = "CLARABEL")
  expect_equal(val, sqrt(2), tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), c(1, 1), tolerance = 1e-5)
})

## @cvxpy test_problem_form.py::TestResolveAndBuildChain::test_solver_not_installed
test_that("ProblemForm parity: unknown solver is rejected", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))

  expect_error(CVXR:::construct_solving_chain(prob, solver = "NONEXISTENT"),
               class = "rlang_error")
})

## @cvxpy test_problem_form.py::TestResolveAndBuildChain::test_solver_cannot_handle
test_that("ProblemForm parity: incompatible solver is rejected", {
  skip_if_not("OSQP" %in% installed_solvers())

  x <- Variable(2)
  prob <- Problem(Minimize(p_norm(x, 2)), list(x >= 1))

  expect_error(CVXR:::construct_solving_chain(prob, solver = "OSQP"),
               "does not support|cannot handle|required cone")
})

## @cvxpy test_problem_form.py::TestGP::test_gp_cones
test_that("ProblemForm parity: GP chain inserts Dgp2Dcp and preserves conic routing", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x + y), list(x * y >= 1))
  chain <- CVXR:::construct_solving_chain(prob, gp = TRUE)

  expect_true(any(vapply(chain@reductions, S7::S7_inherits, logical(1L), Dgp2Dcp)))
  expect_false(any(vapply(chain@reductions, problem_form_is_quad_stuffing, logical(1L))))

  x2 <- Variable(2, pos = TRUE)
  geo_prob <- Problem(Maximize(geo_mean(x2)), list(x2 <= 2))
  geo_chain <- CVXR:::construct_solving_chain(geo_prob, gp = TRUE)
  expect_true(any(vapply(geo_chain@reductions, S7::S7_inherits, logical(1L), Dgp2Dcp)))
})

## @cvxpy test_problem_form.py::TestGP::test_gp_solve_default
test_that("ProblemForm parity: GP with default solver solves correctly", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x + y), list(x * y >= 1))

  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 2, tolerance = 1e-3)
})

## @cvxpy test_problem_form.py::TestGP::test_gp_solve_named
test_that("ProblemForm parity: GP with named SCS solver solves correctly", {
  skip_if_not("SCS" %in% installed_solvers())

  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x + y), list(x * y >= 1))

  expect_equal(problem_form_chain_solver(prob, solver = "SCS", gp = TRUE), "SCS")
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(val, 2, tolerance = 1e-2)
})
