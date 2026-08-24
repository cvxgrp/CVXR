## @cvxpy tests/test_derivative.py
##
## A parameter whose value is 0 must still have a gradient.
##
## Why this file exists: it did not. `gradient(p)` came back as exactly 0
## whenever `value(p) == 0` -- no error, no warning, and correct as soon as the
## parameter moved off zero. Gradient descent started at p = 0 never moves.
##
## Cause, in two halves:
##
##  1. CVXR built `A` by contracting the DPP tensor with the parameter vector.
##     At p = 0 the entry of `A` that p drives multiplies out to zero and DROPS
##     OUT OF THE SPARSITY PATTERN, so there is no entry to differentiate with
##     respect to. CVXPY prevents this with `keep_zeros`
##     (cone_matrix_stuffing.py:213 -> ReducedMat.cache, utilities.py:138-140 ->
##     A_mapping_nonzero_rows), whose single caller is DIFFCP.apply
##     (diffcp_conif.py:72). CVXR had no such flag anywhere.
##
##  2. Even once CVXR emitted the explicit zero, R `diffcp` <= 0.1.1 called
##     `Matrix::drop0(A)` BEFORE reading the pattern `dA` is reported on, so it
##     threw the entry away again. Fixed in diffcp 0.1.2; Python diffcp captures
##     the pattern first (cone_program.py:643-653) and always had.
##
## Both halves are needed, which is why DESCRIPTION requires diffcp >= 0.1.2.

## Closed form: maximize x1 + 2*x2 subject to p*x1 + x2 <= 1, 0 <= x1 <= 3,
## 0 <= x2 <= 5 gives x1 = 3, x2 = 1 - 3p for p in [0, 1/3].
## So  d/dp sum(x)      = -3   and   d/dp (x1 + 2*x2) = -6,  at EVERY p.
.kz_problem <- function(pv) {
  x <- Variable(2)
  p <- Parameter()
  value(p) <- pv
  list(prob = Problem(Maximize(x[1] + 2 * x[2]),
                      list(p * x[1] + x[2] <= 1, x >= 0, x[1] <= 3, x[2] <= 5)),
       x = x, p = p)
}

.kz_grad <- function(pv, seed = NULL) {
  h <- .kz_problem(pv)
  psolve(h$prob, solver = "DIFFCP", requires_grad = TRUE)
  if (!is.null(seed)) gradient(h$x) <- seed
  backward(h$prob)
  as.numeric(gradient(h$p))
}

test_that("a parameter at zero still has a gradient", {
  skip_if_not_installed("diffcp", minimum_version = "0.1.2")
  ## The regression. Pre-fix this returned exactly 0.
  expect_equal(.kz_grad(0), -3, tolerance = 1e-5)
  ## CVXPY 1.9.2 on the same problem returns -3.000000 at p = 0.
})

test_that("the gradient is the same at zero as either side of it", {
  skip_if_not_installed("diffcp", minimum_version = "0.1.2")
  ## The closed form is linear in p, so the gradient is constant. A zero at
  ## p = 0 sitting between two correct values is the exact signature of the bug.
  expect_equal(.kz_grad(0.0), -3, tolerance = 1e-5)
  expect_equal(.kz_grad(0.1), -3, tolerance = 1e-5)
  expect_equal(.kz_grad(0.2), -3, tolerance = 1e-5)
})

test_that("an explicit objective seed gives the objective's gradient", {
  skip_if_not_installed("diffcp", minimum_version = "0.1.2")
  expect_equal(.kz_grad(0.0, seed = c(1, 2)), -6, tolerance = 1e-5)
  expect_equal(.kz_grad(0.1, seed = c(1, 2)), -6, tolerance = 1e-5)
})

test_that("the solve itself is unaffected by keeping zeros", {
  skip_if_not_installed("diffcp")
  h <- .kz_problem(0)
  val <- psolve(h$prob, solver = "DIFFCP", requires_grad = TRUE)
  expect_equal(as.numeric(val), 5, tolerance = 1e-5)
  expect_equal(as.numeric(value(h$x)), c(3, 1), tolerance = 1e-5)
  expect_equal(status(h$prob), "optimal")
})

## The CVXR half on its own, with no solver involved: the data handed to DIFFCP
## must carry the parameter-driven entry as a stored zero. This keeps failing
## loudly if `keep_zeros` is ever dropped from apply_parameters(), even if a
## future diffcp were to paper over it.
test_that("problem_data() for DIFFCP keeps parameter-driven zeros in A", {
  skip_if_not_installed("diffcp")
  h <- .kz_problem(0)
  pd <- problem_data(h$prob, solver = "DIFFCP")
  A <- pd$data[["A"]]
  expect_true(inherits(A, "sparseMatrix"))
  ## At p = 0 that entry is numerically zero but must still be STORED.
  expect_gt(sum(A@x == 0), 0)
  ## ...and it disappears from an ordinary (non-DIFFCP) compile, which is the
  ## contrast that shows the flag is doing the work rather than some accident.
  h2 <- .kz_problem(0)
  pd2 <- problem_data(h2$prob, solver = "SCS")
  expect_equal(sum(pd2$data[["A"]]@x == 0), 0)
})

test_that("A_mapping_nonzero_rows finds the parameter-driven rows", {
  ## CVXPY SOURCE: canonInterface.py:170-180. Unit-level, no solve.
  h <- .kz_problem(0)
  pd <- problem_data(h$prob, solver = "DIFFCP")
  pp <- pd$data[["param_prob"]]
  skip_if(is.null(pp), "solver data carries no param_prob")
  rows <- CVXR:::A_mapping_nonzero_rows(pp@A_tensor, pp@x_length,
                                        const_col = pp@param_id_to_col[["-1"]])
  expect_gt(length(rows), 0)
  expect_true(all(rows >= 1))
  ## Every returned row must lie inside the A block, never the trailing b block.
  expect_true(all(rows <= (nrow(pp@A_tensor) %/% (pp@x_length + 1L)) * pp@x_length))
})
