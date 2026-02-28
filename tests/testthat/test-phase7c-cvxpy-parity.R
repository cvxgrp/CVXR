## Phase 7c CVXPY Parity Tests
## Sources: test_atoms.py, test_domain.py, test_nonlinear_atoms.py, test_constant_atoms.py

# ═══════════════════════════════════════════════════════════════════
# test_atoms.py parity: DCP properties on specific inputs
# ═══════════════════════════════════════════════════════════════════

## CVXPY: test_atoms.py::test_log1p (lines 889-897)
## @cvxpy test_atoms.py::TestAtoms::test_log1p
test_that("CVXPY parity: log1p(1) is NONNEG CONSTANT", {
  expr <- log1p_atom(1)
  expect_true(is_nonneg(expr))
  expect_false(is_nonpos(expr))
  expect_true(is_constant(expr))
  expect_equal(expr@shape, c(1L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_log1p
test_that("CVXPY parity: log1p(-0.5) is NONPOS", {
  expr <- log1p_atom(-0.5)
  expect_false(is_nonneg(expr))
  expect_true(is_nonpos(expr))
})

## CVXPY: test_atoms.py::test_xexp (lines 179-191)
## @cvxpy test_atoms.py::TestAtoms::test_xexp
test_that("CVXPY parity: xexp positive arg → CONVEX NONNEG", {
  x <- Variable(pos = TRUE)
  atom <- xexp(x)
  expect_true(is_convex(atom))
  expect_true(is_nonneg(atom))
})

## @cvxpy test_atoms.py::TestAtoms::test_xexp
test_that("CVXPY parity: xexp negative arg → NOT CONCAVE, NONPOS", {
  x <- Variable(neg = TRUE)
  atom <- xexp(x)
  expect_false(is_concave(atom))
  expect_true(is_nonpos(atom))
})

## CVXPY: test_atoms.py::test_log_sum_exp (lines 1941-1954)
## @cvxpy test_atoms.py::TestAtoms::test_log_sum_exp
test_that("CVXPY parity: log_sum_exp nonneg arg → CONVEX NONNEG", {
  x <- Variable(nonneg = TRUE)
  atom <- log_sum_exp(x)
  expect_true(is_convex(atom))
  expect_true(is_nonneg(atom))
})

## @cvxpy test_atoms.py::TestAtoms::test_log_sum_exp
test_that("CVXPY parity: log_sum_exp nonpos arg → CONVEX UNKNOWN sign", {
  x <- Variable(nonpos = TRUE)
  atom <- log_sum_exp(x)
  expect_true(is_convex(atom))
  ## CVXPY: sign is UNKNOWN (not nonneg, not nonpos)
  expect_false(is_nonneg(atom))
})

# ═══════════════════════════════════════════════════════════════════
# test_constant_atoms.py parity: numeric evaluation
# ═══════════════════════════════════════════════════════════════════

## CVXPY: kl_div(e, 1) = 1, kl_div(e, e) = 0
## @cvxpy NONE
test_that("CVXPY parity: kl_div constant evaluation", {
  e <- exp(1)
  expect_equal(as.numeric(value(kl_div(Constant(e), Constant(1)))),
               1.0, tolerance = 1e-6)
  expect_equal(as.numeric(value(kl_div(Constant(e), Constant(e)))),
               0.0, tolerance = 1e-6)
  ## Vector: kl_div([e, 1], 1) = [1, 0]
  result <- as.numeric(value(kl_div(Constant(c(e, 1)), Constant(1))))
  expect_equal(result, c(1, 0), tolerance = 1e-6)
})

## CVXPY: rel_entr(e, 1) = e, rel_entr(e, e) = 0
## @cvxpy NONE
test_that("CVXPY parity: rel_entr constant evaluation", {
  e <- exp(1)
  expect_equal(as.numeric(value(rel_entr(Constant(e), Constant(1)))),
               e, tolerance = 1e-6)
  expect_equal(as.numeric(value(rel_entr(Constant(e), Constant(e)))),
               0.0, tolerance = 1e-6)
  ## Vector: rel_entr([e, 1], 1) = [e, 0]
  result <- as.numeric(value(rel_entr(Constant(c(e, 1)), Constant(1))))
  expect_equal(result, c(e, 0), tolerance = 1e-6)
})

## CVXPY: logistic(log(5)) = log(6), logistic(0) = log(2), etc.
## @cvxpy NONE
test_that("CVXPY parity: logistic constant evaluation", {
  ## logistic(x) = log(1 + exp(x))
  logistic_in <- matrix(c(log(5), 0, log(7), log(0.3)), nrow = 2, ncol = 2)
  logistic_expected <- matrix(c(log(6), log(2), log(8), log(1.3)), nrow = 2, ncol = 2)
  result <- as.matrix(value(logistic(Constant(logistic_in))))
  expect_equal(result, logistic_expected, tolerance = 1e-6)
})

## CVXPY: log_sum_exp([[5,7],[0,-3]]) ≈ 7.1277708269
## @cvxpy NONE
test_that("CVXPY parity: log_sum_exp constant evaluation", {
  ## Note: R column-major so matrix(c(5,0,7,-3), 2, 2) = [[5,7],[0,-3]] in row-major
  m <- matrix(c(5, 0, 7, -3), nrow = 2, ncol = 2)
  result <- as.numeric(value(log_sum_exp(Constant(m))))
  expect_equal(result, 7.1277708269, tolerance = 1e-6)
})

## CVXPY: log_sum_exp axis=0 on [[5,7,1],[0,-3,6]]
## @cvxpy NONE
test_that("CVXPY parity: log_sum_exp axis=0 constant evaluation", {
  ## R: matrix(c(5,0,7,-3,1,6), 2, 3) = row0=[5,7,1], row1=[0,-3,6]
  m <- matrix(c(5, 0, 7, -3, 1, 6), nrow = 2, ncol = 3)
  result <- as.numeric(value(log_sum_exp(Constant(m), axis = 2L, keepdims = TRUE)))
  ## CVXPY keepdims=True, axis=0 gives [[5.00671535, 7.0000454, 6.00671535]]
  expect_equal(result, c(5.00671535, 7.0000454, 6.00671535), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("CVXPY parity: log_sum_exp axis=1 constant evaluation", {
  m <- matrix(c(5, 0, 7, -3, 1, 6), nrow = 2, ncol = 3)
  result <- as.numeric(value(log_sum_exp(Constant(m), axis = 1L)))
  ## CVXPY axis=1 gives [7.12910891, 6.00259878]
  expect_equal(result, c(7.12910891, 6.00259878), tolerance = 1e-4)
})

## CVXPY: tv([1,-1,2]) = 5
## @cvxpy NONE
test_that("CVXPY parity: total_variation constant 1D", {
  expect_equal(as.numeric(value(total_variation(Constant(c(1, -1, 2))))),
               5.0, tolerance = 1e-6)
})

## CVXPY: tv([[-5,2],[-3,1]]) = sqrt(53)
## @cvxpy NONE
test_that("CVXPY parity: total_variation constant 2D", {
  m <- matrix(c(-5, -3, 2, 1), nrow = 2, ncol = 2)
  expect_equal(as.numeric(value(total_variation(Constant(m)))),
               sqrt(53), tolerance = 1e-4)
})

## CVXPY: tv([[3,4,5],[6,7,8],[9,10,11]]) = 4*sqrt(10)
## @cvxpy NONE
test_that("CVXPY parity: total_variation constant 2D uniform gradient", {
  m <- matrix(c(3, 6, 9, 4, 7, 10, 5, 8, 11), nrow = 3, ncol = 3)
  expect_equal(as.numeric(value(total_variation(Constant(m)))),
               4 * sqrt(10), tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════
# test_domain.py parity: domain constraint solving
# ═══════════════════════════════════════════════════════════════════

## CVXPY: test_domain.py::test_log1p — minimize a s.t. log1p(a) domain → a = -1
## @cvxpy test_domain.py::TestDomain::test_log1p
test_that("CVXPY parity: log1p domain boundary", {
  skip_if_not_installed("clarabel")
  a <- Variable()
  dom <- domain(log1p_atom(a))
  prob <- Problem(Minimize(a), dom)
  psolve(prob, solver = "CLARABEL", verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(a)), -1, tolerance = 1e-4)
})

## CVXPY: test_domain.py::test_kl_div — minimize a+b s.t. kl_div(a,b) domain → a=0, b=0
## @cvxpy test_domain.py::TestDomain::test_kl_div
test_that("CVXPY parity: kl_div domain boundary", {
  skip_if_not_installed("clarabel")
  a <- Variable()
  b <- Variable()
  dom <- domain(kl_div(a, b))
  prob <- Problem(Minimize(a + b), dom)
  psolve(prob, solver = "CLARABEL", verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(a)), 0, tolerance = 1e-3)
  expect_equal(as.numeric(value(b)), 0, tolerance = 1e-3)
})

## CVXPY: test_domain.py::test_rel_entr — minimize a+b s.t. rel_entr(a,b) domain → a=0, b=0
## @cvxpy test_domain.py::TestDomain::test_rel_entr
test_that("CVXPY parity: rel_entr domain boundary", {
  skip_if_not_installed("clarabel")
  a <- Variable()
  b <- Variable()
  dom <- domain(rel_entr(a, b))
  prob <- Problem(Minimize(a + b), dom)
  psolve(prob, solver = "CLARABEL", verbose = FALSE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(a)), 0, tolerance = 1e-3)
  expect_equal(as.numeric(value(b)), 0, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════
# test_nonlinear_atoms.py parity: kl_div vs rel_entr difference
# ═══════════════════════════════════════════════════════════════════

## CVXPY: test_nonlinear_atoms.py::test_difference_kl_div_rel_entr (lines 123-143)
## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_difference_kl_div_rel_entr
test_that("CVXPY parity: kl_div vs rel_entr difference", {
  skip_if_not_installed("clarabel")

  x <- Variable()
  y <- Variable()

  ## kl_div minimization: x == y at optimum
  kl_prob <- Problem(Minimize(kl_div(x, y)), list(x + y <= 1))
  kl_val <- psolve(kl_prob, solver = "CLARABEL", verbose = FALSE)
  expect_equal(status(kl_prob), "optimal")
  ## kl_div(x,y) = 0 when x == y
  expect_equal(kl_val, 0.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), as.numeric(value(y)), tolerance = 1e-3)

  ## rel_entr minimization: asymmetric solution
  ## Wolfram Alpha: x ≈ 0.2178, y ≈ 0.7822, obj ≈ -0.2785
  rel_prob <- Problem(Minimize(rel_entr(x, y)), list(x + y <= 1))
  rel_val <- psolve(rel_prob, solver = "CLARABEL", verbose = FALSE)
  expect_equal(status(rel_prob), "optimal")
  expect_equal(as.numeric(value(x)), 0.2178117, tolerance = 1e-3)
  expect_equal(as.numeric(value(y)), 0.7821882, tolerance = 1e-3)
  expect_equal(rel_val, -0.278464, tolerance = 1e-3)
})
