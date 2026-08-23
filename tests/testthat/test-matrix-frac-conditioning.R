## @cvxpy tests/test_atoms.py::TestAtoms::test_matrix_frac
##
## `matrix_frac()` must return a value for an ILL-CONDITIONED second argument,
## as CVXPY does, and must say so.
##
## Every inverse in matrix_frac.py -- numeric (:42,44), _grad (:68) and the
## constructor shortcut (:149) -- is `numpy.linalg.inv`, which has no
## conditioning test and raises only on EXACT singularity. Base R's `solve()`
## applies `tol = .Machine$double.eps` against the reciprocal condition number
## and refuses below it. So for a merely ill-conditioned P the two languages
## disagreed about whether a value exists at all.
##
## Found by revdepcheck, not by this suite: SLSEdesign's vignette computes a
## c-optimal design whose optimum is DEGENERATE -- the information matrix is
## rank-deficient in the directions the criterion does not constrain -- so the
## better the solver converges, the closer P gets to singular. At the optimum
## rcond(P) = 4.3e-18, `solve(P)` refused, and that aborted
## `value(problem@objective)` inside psolve(), taking the whole solve with it.
## CVXR 1.9.1 did not hit it only because it stopped further out, where it
## returned -4.08e-05 -- also wrong, since matrix_frac is mathematically >= 0,
## and four orders worse than the -5.2e-14 the converged answer gives.
##
## CVXR deviates from upstream in ONE way, deliberately: it warns. numpy
## inverts silently and the result can be meaningless; R's refusal existed for
## a reason. Keeping the value and adding the warning gives both.

## The exact matrix from that vignette's optimum, to 17 significant digits.
.ill_conditioned_B <- function() {
  matrix(c(1, 0, 0,
           0, 3.8709486609044189e+12, 2.2715196466934160e+15,
           0, 2.2715196466934160e+15, 2.2853429570577206e+17),
         nrow = 3, byrow = TRUE)
}

test_that("the fixture really is below base solve()'s tolerance", {
  B <- .ill_conditioned_B()
  expect_lt(rcond(B), .Machine$double.eps)
  ## base solve() refuses it -- this is the behavior being deviated from
  expect_error(solve(B), "computationally singular")
})

test_that("matrix_frac evaluates an ill-conditioned P, and warns", {
  B <- .ill_conditioned_B()
  expr <- matrix_frac(Constant(matrix(c(0, 1, 1), ncol = 1)), Constant(B))
  expect_warning(value(expr), class = "matrixFracConditioning")
  val <- suppressWarnings(value(expr))
  ## numpy on the identical matrix: -5.2394597073595286e-14
  expect_equal(as.numeric(val), -5.2394597073595286e-14, tolerance = 1e-6)
})

test_that("a well-conditioned P neither warns nor changes answer", {
  P <- diag(c(1, 2, 4))
  X <- Constant(matrix(c(1, 1, 1), ncol = 1))
  expr <- matrix_frac(X, Constant(P))
  expect_no_warning(v <- value(expr))
  ## t(x) P^-1 x = 1 + 1/2 + 1/4
  expect_equal(as.numeric(v), 1.75, tolerance = 1e-12)
})

test_that("an exactly singular P is still an error, with no extra warning", {
  ## LA.inv raises LinAlgError here too; upstream never decided what the value
  ## should be (matrix_frac.py:38 is a standing TODO), so CVXR keeps erroring.
  ## And it must not ALSO warn "ill-conditioned" on the way there -- one problem,
  ## one message.
  P <- matrix(c(1, 0, 0, 0), 2, 2)
  X <- Constant(matrix(c(1, 1), ncol = 1))
  expect_error(value(matrix_frac(X, Constant(P))), "singular")
  expect_no_warning(tryCatch(value(matrix_frac(X, Constant(P))),
                             error = function(e) NULL))
})

test_that("the constructor shortcut for a numeric P behaves the same way", {
  ## CVXPY SOURCE: matrix_frac.py:147-152 -- a plain numeric P is inverted
  ## eagerly into a quad_form.
  x <- Variable(3)
  expect_no_error(matrix_frac(x, diag(c(1, 2, 4))))
  expect_no_warning(matrix_frac(x, diag(c(1, 2, 4))))
  ## exactly singular -> error, and no conditioning warning alongside it
  expect_error(matrix_frac(x, matrix(0, 3, 3)), "singular")
  expect_no_warning(tryCatch(matrix_frac(x, matrix(0, 3, 3)),
                             error = function(e) NULL))
})

test_that("the SLSEdesign shape solves end to end", {
  ## The reduced form of SLSEdesign::copt()'s problem: a c-optimal design whose
  ## optimum drives the information matrix to rank deficiency. Pre-fix this
  ## raised "system is computationally singular" from inside psolve().
  my_peleg <- function(xi, theta) {
    deno <- theta[1] + theta[2] * xi
    matrix(c(-xi / deno^2, -xi^2 / deno^2), ncol = 1)
  }
  N <- 201
  u <- seq(0, 100, length.out = N)
  theta <- c(0.5, 0.05)
  w <- Variable(N)
  mf <- sapply(u, my_peleg, theta)
  G2 <- mf %*% DiagVec(w) %*% t(mf)
  B <- bmat(list(list(1, matrix(0, 1, 2)), list(matrix(0, 2, 1), G2)))
  prob <- Problem(Minimize(matrix_frac(c(0, 1, 1), B)),
                  list(w >= 0, sum_entries(w) == 1))
  val <- suppressWarnings(psolve(prob, solver = "SCS", reltol = 1e-6, abstol = 1e-6))
  expect_true(is.numeric(val))
  expect_false(is.na(as.numeric(val)))
  ## The true optimum is 0; SCS approaches it from below on this degenerate
  ## problem. Pin only that it is near zero, not its sign.
  expect_lt(abs(as.numeric(val)), 1e-3)
})
