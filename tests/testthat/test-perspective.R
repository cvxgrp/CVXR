## Tests for Perspective atom and canonicalizer
## CVXPY reference: cvxpy/tests/test_perspective.py

## @cvxpy test_perspective.py::test_monotonicity
test_that("Perspective properties: convex f gives convex perspective", {
  x <- Variable(nonneg = TRUE)
  s <- Variable(nonneg = TRUE)
  f <- square(x)
  p <- Perspective(f, s)

  expect_true(is_convex(p))
  expect_false(is_concave(p))
  expect_equal(p@shape, c(1L, 1L))
})

## @cvxpy test_perspective.py::test_monotonicity
test_that("Perspective properties: concave f gives concave perspective", {
  x <- Variable(nonneg = TRUE)
  s <- Variable(nonneg = TRUE)
  f <- log(x)
  p <- Perspective(f, s)

  expect_false(is_convex(p))
  expect_true(is_concave(p))
  expect_equal(p@shape, c(1L, 1L))
})

## @cvxpy test_perspective.py::test_scalar_x
test_that("Perspective sign: nonneg f and nonneg s", {
  x <- Variable(nonneg = TRUE)
  s <- Variable(nonneg = TRUE)
  f <- square(x)  # nonneg
  p <- Perspective(f, s)

  expect_true(is_nonneg(p))
  expect_false(is_nonpos(p))
})

## @cvxpy test_perspective.py::test_afine_s
test_that("Perspective is_dcp works", {
  x <- Variable(nonneg = TRUE)
  s <- Variable(nonneg = TRUE)

  ## Convex f -> Perspective is convex -> DCP
  f_conv <- square(x)
  p_conv <- Perspective(f_conv, s)
  expect_true(is_dcp(Minimize(p_conv)))

  ## Concave f -> Perspective is concave -> DCP in Maximize
  f_conc <- log(x)
  p_conc <- Perspective(f_conc, s)
  expect_true(is_dcp(Maximize(p_conc)))
})

## @cvxpy test_perspective.py::test_evaluate_persp
test_that("Perspective numeric evaluation: s * f(x/s)", {
  x <- Variable(nonneg = TRUE)
  s <- Variable(nonneg = TRUE)
  f <- square(x)
  p <- Perspective(f, s)

  value(x) <- 3.0
  value(s) <- 2.0
  ## s * f(x/s) = 2 * (3/2)^2 = 2 * 2.25 = 4.5
  expect_equal(value(p), 4.5, tolerance = 1e-8)
})

## @cvxpy test_perspective.py::test_evaluate_persp
test_that("Perspective numeric evaluation: concave log", {
  x <- Variable(nonneg = TRUE)
  s <- Variable(nonneg = TRUE)
  f <- log(x)
  p <- Perspective(f, s)

  value(x) <- 4.0
  value(s) <- 2.0
  ## s * log(x/s) = 2 * log(4/2) = 2 * log(2) = 1.3862944
  expect_equal(value(p), 2 * log(2), tolerance = 1e-8)
})

## @cvxpy test_perspective.py::test_assert_s_nonzero
test_that("Perspective validation: non-scalar f errors", {
  x <- Variable(3, nonneg = TRUE)
  s <- Variable(nonneg = TRUE)
  expect_error(Perspective(x, s), "scalar")
})

## @cvxpy test_perspective.py::test_assert_s_nonzero
test_that("Perspective validation: s not nonneg errors", {
  x <- Variable(nonneg = TRUE)
  s <- Variable()  # not nonneg
  f <- square(x)
  expect_error(Perspective(f, s), "nonneg")
})

## @cvxpy test_perspective.py::test_assert_s_nonzero
test_that("Perspective validation: s not Variable errors", {
  x <- Variable(nonneg = TRUE)
  f <- square(x)
  expect_error(Perspective(f, Constant(1)), "Variable")
})

## @cvxpy test_perspective.py::test_power
test_that("Perspective solve: minimize convex square", {
  skip_if_not_installed("clarabel")
  x <- Variable(nonneg = TRUE)
  s <- Variable(nonneg = TRUE)
  f <- square(x)
  obj <- Minimize(perspective(f, s))
  prob <- Problem(obj, list(x + s == 1, s >= 0.1))
  opt_val <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), "optimal")
  ## Optimal: x -> 0, s -> 1, value -> 0
  expect_equal(opt_val, 0, tolerance = 1e-3)
})

## @cvxpy test_perspective.py::test_exp
test_that("Perspective solve: maximize concave log", {
  skip_if_not_installed("clarabel")
  x <- Variable(nonneg = TRUE)
  s <- Variable(nonneg = TRUE)
  f <- log(x)
  obj <- Maximize(perspective(f, s))
  prob <- Problem(obj, list(x + s <= 2, x >= 0.5, s >= 0.1))
  opt_val <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), "optimal")
  ## CVXPY reference: 0.5569290849
  expect_equal(opt_val, 0.556929, tolerance = 1e-3)
})

## @cvxpy test_perspective.py::test_exp
test_that("Perspective solve: minimize exp", {
  skip_if_not_installed("clarabel")
  x <- Variable()
  s <- Variable(nonneg = TRUE)
  f <- exp(x)
  obj <- Minimize(perspective(f, s))
  prob <- Problem(obj, list(x >= 1, s >= 0.5, s <= 2))
  opt_val <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), "optimal")
  ## Optimal at x=1, s=1: s*exp(x/s) = exp(1) = 2.71828...
  expect_equal(opt_val, exp(1), tolerance = 1e-2)
})

## @cvxpy test_perspective.py::test_quad_atom
test_that("Perspective solve: quad_form with bounded s", {
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  s <- Variable(nonneg = TRUE)
  P <- matrix(c(2, 1, 1, 3), 2, 2)
  f <- quad_form(x, P)
  obj <- Minimize(perspective(f, s))
  prob <- Problem(obj, list(x[1] == 1, x[2] == 1, s >= 0.5, s <= 2))
  opt_val <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), "optimal")
  ## f(x) = x'Px = 7, perspective = 7/s, minimized at s=2 -> 3.5
  expect_equal(opt_val, 3.5, tolerance = 1e-2)
})

## @cvxpy test_perspective.py::test_rel_entr
test_that("Perspective solve: rel_entr multi-variable f", {
  skip_if_not_installed("clarabel")
  x <- Variable(nonneg = TRUE)
  y <- Variable(nonneg = TRUE)
  s <- Variable(nonneg = TRUE)
  f <- rel_entr(x, y)
  obj <- Minimize(perspective(f, s))
  prob <- Problem(obj, list(x == 2, y == 1, s >= 0.5, s <= 3))
  opt_val <- psolve(prob, solver = "CLARABEL")

  expect_equal(status(prob), "optimal")
  ## perspective(rel_entr(x,y), s) = x*log(x/y) = 2*log(2) regardless of s
  expect_equal(opt_val, 2 * log(2), tolerance = 1e-3)
})

## @cvxpy test_perspective.py::test_power
test_that("Perspective solve: SCS solver", {
  skip_if_not_installed("scs")
  x <- Variable(nonneg = TRUE)
  s <- Variable(nonneg = TRUE)
  f <- square(x)
  obj <- Minimize(perspective(f, s))
  prob <- Problem(obj, list(x + s == 1, s >= 0.1))
  opt_val <- psolve(prob, solver = "SCS")

  expect_equal(status(prob), "optimal")
  expect_equal(opt_val, 0, tolerance = 1e-2)
})

## @cvxpy test_perspective.py::test_scalar_x
test_that("Perspective: user-facing perspective() function", {
  x <- Variable(nonneg = TRUE)
  s <- Variable(nonneg = TRUE)
  f <- square(x)
  p <- perspective(f, s)

  expect_true(S7_inherits(p, Perspective))
  expect_true(is_convex(p))
})

## @cvxpy test_perspective.py::test_monotonicity
test_that("Perspective: monotonicity is FALSE", {
  x <- Variable(nonneg = TRUE)
  s <- Variable(nonneg = TRUE)
  f <- square(x)
  p <- Perspective(f, s)

  expect_false(is_incr(p, 1L))
  expect_false(is_decr(p, 1L))
})

## perspective_canon links a symmetric/PSD variable to its reduced upper
## triangle. CVXPY pairs `var[np.triu_indices(n)]` (ROW-major,
## perspective_canon.py:97-100), the order CvxAttr2Constr gives the reduced
## variable through upper_tri_to_full. CVXR paired the COLUMN-major triangle
## until 1.9.2.9140, which coincides for n = 2 and swaps the (1,3)/(2,2)
## entries at n = 3 (claims-audit row C0642, repro 23). CVXPY's own PSD tests
## use targets (identity, trace) that this swap leaves invariant, so the
## distinguishing target here is R-specific: upper triangle numbered 1..6.
## @cvxpy NONE
test_that("Perspective: symmetric and PSD variables pair the row-major upper triangle", {
  skip_if_not_installed("clarabel")
  C <- matrix(c(1, 2, 3, 2, 4, 5, 3, 5, 6), 3, 3)  # C[1, 3] = 3 != C[2, 2] = 4
  s_fun <- function(X, target) {
    s <- Variable(nonneg = TRUE)
    prob <- Problem(Minimize(perspective(sum_squares(X - target), s)), list(s == 1))
    val <- psolve(prob, solver = "CLARABEL")
    expect_equal(status(prob), "optimal")
    expect_equal(val, 0, tolerance = 1e-6)
    expect_equal(value(X), target, tolerance = 1e-5)
  }
  s_fun(Variable(c(3, 3), symmetric = TRUE), C)
  s_fun(Variable(c(3, 3), PSD = TRUE), C + 9 * diag(3))
})

## The PSD cone comes from an ATOM inside f, which the outer chain's cone scan
## does not see, so the PSD -> SvecPSD conversion has to happen in the
## perspective's own aux chain, built for the outer solver as CVXPY's
## perspective_canon.py:35-38 does, and perspective_utils.py:55-56 rebuilds
## the SvecPSD. CVXR built that chain without the solver until 1.9.2.9177, and
## both solvers here failed on a cone-dimension mismatch. CVXPY's own PSD
## perspective tests use PSD variables, which the outer scan does see.
## CVXPY 1.9.3: 3.618034 (= lambda_max(C)) on CLARABEL and SCS (repro 36).
## @cvxpy NONE
test_that("Perspective of a PSD-cone atom solves on svec-PSD solvers", {
  C <- matrix(c(2, 1, 1, 3), 2, 2)
  expected <- max(eigen(C, symmetric = FALSE, only.values = TRUE)$values)
  for (solver in c("CLARABEL", "SCS")) {
    if (!(solver %in% installed_solvers())) next
    X <- Variable(c(2, 2), symmetric = TRUE)
    s <- Variable(nonneg = TRUE)
    prob <- Problem(Minimize(perspective(lambda_max(X), s)), list(X == C, s == 2))
    val <- psolve(prob, solver = solver)
    expect_equal(status(prob), "optimal")
    expect_equal(val, expected, tolerance = 1e-4)
  }
})

## @cvxpy NONE
test_that("Perspective graph_implementation errors", {
  x <- Variable(nonneg = TRUE)
  s <- Variable(nonneg = TRUE)
  f <- square(x)
  p <- Perspective(f, s)

  expect_error(graph_implementation(p, list(), c(1L, 1L)), "Dcp2Cone")
})
