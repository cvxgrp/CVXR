## CVXPY 1.9.0 PR #3180 parity tests (multi-bug "robustness" PR)
## Ported from: cvxpy/tests (CVXPY 1.9.0, worktree ~/GitHub/cvxpy-1.9)
##
## #3180 is a multi-bug PR; each fix landed a small new test in a different
## CVXPY test file. CVXR already carried several of these behaviors (ported
## from 1.8.2); the rest are small source fixes landed alongside this file:
##   - power() const-base / variable-exponent dispatch (power.R, expression.R)
##   - Power/PowerApprox reject a non-Constant/Parameter exponent (power.R)
##   - diag(v).is_hermitian requires real v (diag.R)
##   - Clarabel "Unsolved" status -> SOLVER_ERROR (clarabel_conif.R)
##
## The two quantum-atom tests (quantum_rel_entr, von_neumann_entr) are skip()ed
## and @v19-pending: those atoms are DEFERRED post-v1.0 (see CLAUDE.md /
## atoms/quantum_rel_entr.R, von_neumann_entr.R stubs).

# =====================================================================
# test_complex.py::TestComplex2RealAccepts::test_accepts_returns_bool
# =====================================================================
## @cvxpy test_complex.py::TestComplex2RealAccepts::test_accepts_returns_bool
test_that("CVXPY 1.9 parity (#3180): Complex2Real accepts() returns a bool", {
  red <- Complex2Real()

  ## Complex problem should be accepted.
  x <- Variable(c(2L, 2L), complex = TRUE)
  prob <- Problem(Minimize(norm(x, "F")), list(x == diag(2)))
  acc <- reduction_accepts(red, prob)
  expect_true(is.logical(acc) && length(acc) == 1L)
  expect_identical(acc, TRUE)

  ## Real problem should not be accepted.
  y <- Variable(c(2L, 2L))
  prob_real <- Problem(Minimize(norm(y, "F")), list(y == diag(2)))
  acc_real <- reduction_accepts(red, prob_real)
  expect_true(is.logical(acc_real) && length(acc_real) == 1L)
  expect_identical(acc_real, FALSE)
})

# =====================================================================
# test_complex.py::TestComplex::test_diag_vec_hermitian
# =====================================================================
## @cvxpy test_complex.py::TestComplex::test_diag_vec_hermitian
test_that("CVXPY 1.9 parity (#3180): diag(v) is Hermitian only when v is real", {
  v_complex <- Variable(3L, complex = TRUE)
  D_complex <- diag(v_complex)
  expect_false(is_hermitian(D_complex))

  v_real <- Variable(3L)
  D_real <- diag(v_real)
  expect_true(is_hermitian(D_real))
})

# =====================================================================
# test_constraints.py::TestConstraints::test_dual_cone_no_args
# =====================================================================
## @cvxpy test_constraints.py::TestConstraints::test_dual_cone_no_args
test_that("CVXPY 1.9 parity (#3180): dual_cone() without args works for PSD and SOC", {
  skip_if_not_installed("clarabel")

  ## PSD
  X <- Variable(c(2L, 2L), symmetric = TRUE)
  psd_constr <- PSD(X)
  prob <- Problem(Minimize(matrix_trace(X)), list(psd_constr, matrix_trace(X) >= 1))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  dual_psd <- dual_cone(psd_constr)
  expect_true(S7::S7_inherits(dual_psd, PSD))

  ## SOC
  xx <- Variable(3L)
  tt <- Variable(1L)
  soc_constr <- SOC(tt, xx)
  prob2 <- Problem(Minimize(tt), list(soc_constr, xx == rep(1, 3)))
  psolve(prob2, solver = "CLARABEL")
  expect_equal(status(prob2), "optimal")
  dual_soc <- dual_cone(soc_constr)
  expect_true(S7::S7_inherits(dual_soc, SOC))
})

# =====================================================================
# test_problem.py::TestProblem::test_canonicalization_invert_none_duals
# =====================================================================
## @cvxpy test_problem.py::TestProblem::test_canonicalization_invert_none_duals
test_that("CVXPY 1.9 parity (#3180): Canonicalization.invert handles absent dual_vars", {
  x <- Variable(2L)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  red <- Dcp2Cone()
  res <- reduction_apply(red, prob)
  inv_data <- res[[2L]]

  ## Simulate a solver returning no dual_vars (e.g. MIP solvers): in R that is
  ## an empty list (CVXPY's None). invert() must not error and must yield no
  ## recovered duals.
  sol <- Solution(status = "optimal", opt_val = 2.0,
                  primal_vars = list(), dual_vars = list(), attr = list())
  inv <- reduction_invert(red, sol, inv_data)
  expect_equal(length(inv@dual_vars), 0L)
})

# =====================================================================
# test_clarabel.py::ClarabelTest::test_unsolved_status_handling
# =====================================================================
## @cvxpy test_clarabel.py::ClarabelTest::test_unsolved_status_handling
test_that("CVXPY 1.9 parity (#3180): Clarabel 'Unsolved' status maps to SOLVER_ERROR", {
  ## Map entry (the actual fix).
  expect_equal(CLARABEL_STATUS_MAP[["Unsolved"]], SOLVER_ERROR)

  ## invert() with an "Unsolved" status must not raise and must report
  ## SOLVER_ERROR. R clarabel returns the status as an integer index into
  ## solver_status_descriptions(); index 1 == "Unsolved".
  unsolved_idx <- match("Unsolved", names(clarabel::solver_status_descriptions()))
  mock_solution <- list(status = unsolved_idx, solve_time = 0.01,
                        iterations = 0L, x = NULL, z = NULL)
  sol <- reduction_invert(Clarabel_Solver(), mock_solution, NULL)
  expect_equal(sol@status, SOLVER_ERROR)
})

# =====================================================================
# test_expressions.py::TestExpressions::test_parameter_infinite_values
# =====================================================================
## @cvxpy test_expressions.py::TestExpressions::test_parameter_infinite_values
test_that("CVXPY 1.9 parity (#3180): Parameter accepts/rejects +-Inf by sign", {
  ## Unconstrained scalar parameter accepts +Inf and -Inf.
  p <- Parameter()
  value(p) <- Inf
  expect_equal(as.numeric(value(p)), Inf)
  value(p) <- -Inf
  expect_equal(as.numeric(value(p)), -Inf)

  ## Nonneg parameter accepts +Inf but rejects -Inf.
  p <- Parameter(nonneg = TRUE)
  value(p) <- Inf
  expect_equal(as.numeric(value(p)), Inf)
  expect_error(value(p) <- -Inf, "nonneg")

  ## Nonpos parameter accepts -Inf but rejects +Inf.
  p <- Parameter(nonpos = TRUE)
  value(p) <- -Inf
  expect_equal(as.numeric(value(p)), -Inf)
  expect_error(value(p) <- Inf, "nonpos")

  ## Vector parameter with mixed finite and infinite values.
  p <- Parameter(3L)
  value(p) <- c(1.0, Inf, -Inf)
  expect_equal(as.numeric(value(p)), c(1.0, Inf, -Inf))
})

# =====================================================================
# test_expressions.py::TestExpressions::test_power_const_base
# =====================================================================
## @cvxpy test_expressions.py::TestExpressions::test_power_const_base
test_that("CVXPY 1.9 parity (#3180): power(b, x) with constant/parametric base b", {
  skip_if_not_installed("clarabel")
  x <- Variable()

  ## Basic expression created correctly.
  expect_false(is.null(power(2, x)))

  ## minimize power(2, x) with x in [1, 3] -> x = 1, value = 2.
  prob <- Problem(Minimize(power(2, x)), list(x >= 1, x <= 3))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(prob)), 2.0, tolerance = 1e-3)

  ## base = 10, x >= 2 -> value = 100.
  prob2 <- Problem(Minimize(power(10, x)), list(x >= 2))
  psolve(prob2, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(prob2)), 100.0, tolerance = 1e-2)

  ## Negative base is rejected.
  expect_error(power(-2, x), "positive")

  ## Parameter(pos = TRUE) as base.
  b <- Parameter(pos = TRUE)
  value(b) <- 2.0
  expect_false(is.null(power(b, x)))
  prob3 <- Problem(Minimize(power(b, x)), list(x >= 1, x <= 3))
  psolve(prob3, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-3)
})

# =====================================================================
# test_expressions.py::TestExpressions::test_rpow
# =====================================================================
## @cvxpy test_expressions.py::TestExpressions::test_rpow
test_that("CVXPY 1.9 parity (#3180): a^x (rpow) for constant base a", {
  skip_if_not_installed("clarabel")
  x <- Variable()

  ## Basic expression.
  expect_false(is.null(2^x))

  ## minimize 2^x with x in [1, 3] -> x = 1, value = 2.
  prob <- Problem(Minimize(2^x), list(x >= 1, x <= 3))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(prob)), 2.0, tolerance = 1e-3)

  ## base = 10, x >= 2 -> value = 100.
  prob2 <- Problem(Minimize(10^x), list(x >= 2))
  psolve(prob2, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(prob2)), 100.0, tolerance = 1e-2)

  ## base = 1 gives value 1 regardless of x.
  prob3 <- Problem(Minimize(1^x), list(x >= 0))
  psolve(prob3, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob3)), 1.0, tolerance = 1e-3)

  ## Negative base is rejected.
  expect_error((-2)^x, "positive")

  ## Non-integer base 0.5^x, minimize over x in [-2, 0] -> x = 0, value = 1.
  prob4 <- Problem(Minimize(0.5^x), list(x >= -2, x <= 0))
  psolve(prob4, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 0.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(prob4)), 1.0, tolerance = 1e-3)

  ## Non-constant base is rejected (y^x with y a Variable).
  y <- Variable()
  expect_error(y^x, "Constant or a Parameter")
})

# =====================================================================
# test_quantum_rel_entr.py::TestQuantumRelEntr::test_real_inputs
# =====================================================================
## @cvxpy test_quantum_rel_entr.py::TestQuantumRelEntr::test_real_inputs
## @v19-pending: quantum atoms (quantum_rel_entr) [deferred-post-v1.0]
test_that("CVXPY 1.9 parity (#3180): quantum_rel_entr with real inputs", {
  skip("@v19-pending: quantum atoms deferred to post-v1.0 (see atoms/quantum_rel_entr.R stub)")
  n <- 2L
  X <- Variable(c(n, n), symmetric = TRUE)
  Y <- diag(n)
  prob <- Problem(Minimize(quantum_rel_entr(X, Y)),
                  list(matrix_trace(X) == 1, PSD(X)))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(prob)), -log(n), tolerance = 1e-3)
})

# =====================================================================
# test_von_neumann_entr.py::Test_von_neumann_entr::test_max_entropy
# =====================================================================
## @cvxpy test_von_neumann_entr.py::Test_von_neumann_entr::test_max_entropy
## @v19-pending: quantum atoms (von_neumann_entr) [deferred-post-v1.0]
test_that("CVXPY 1.9 parity (#3180): von_neumann_entr max entropy is log(n)", {
  skip("@v19-pending: quantum atoms deferred to post-v1.0 (see atoms/von_neumann_entr.R stub)")
  n <- 3L
  N <- Variable(c(n, n), symmetric = TRUE)
  prob <- Problem(Maximize(von_neumann_entr(N)),
                  list(matrix_trace(N) == 1, PSD(N)))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(prob)), log(n), tolerance = 1e-3)
})
