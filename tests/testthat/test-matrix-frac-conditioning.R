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
## Found by revdepcheck, not by this suite: SLSEdesign's vignette solves a
## c-optimal design with SCS, and SCS NEVER converges on it (100,000
## iterations, `optimal_inaccurate`, on every platform). Its unconverged
## iterate wanders through near-singular territory: at SCS's stall point
## rcond(B) = 4.3e-18, `solve(B)` refused, and that aborted
## `value(problem@objective)` inside psolve(), taking the whole solve with it.
## The OPTIMUM itself is well-conditioned: CLARABEL and MOSEK independently
## agree on 0.016494 (`optimal`), where rcond(B) = 5.9e-05 and the 2x2
## information block has eigenvalues {16481, 56} -- full rank, two support
## points for two parameters (measured 2026-08-23). An earlier version of
## this header claimed the optimum was degenerate and ~0; both claims were
## derived from unconverged SCS output (CVXR 1.9.1's -4.08e-05 was the same
## artifact) and were false.
##
## CVXR deviates from upstream in ONE way, deliberately: it warns. numpy
## inverts silently and the result can be meaningless; R's refusal existed for
## a reason. Keeping the value and adding the warning gives both.

## The exact matrix at SCS's unconverged stall point on that vignette's
## problem, to 17 significant digits. NOT the optimum -- see the header.
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
  ## The reduced form of SLSEdesign::copt()'s problem. Pre-fix, SCS's
  ## unconverged iterate made B ill-conditioned enough that psolve() raised
  ## "system is computationally singular" from inside value().
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
  ## The SCS guard is ONLY that the solve completes and returns a number.
  ## SCS never converges here (100k iterations, `optimal_inaccurate` on
  ## every platform) and its terminal iterate is BLAS-dependent (-3.3e-09
  ## under Apple Accelerate, 0.096 under reference BLAS), so nothing about
  ## its value may be asserted. An earlier version pinned |val| < 1e-3 on
  ## the false premise that the optimum is 0 -- see the header.
  expect_true(is.numeric(val))
  expect_false(is.na(as.numeric(val)))

  ## The real pin, on a solver that terminates: the true optimum is
  ## 0.0164945 (CLARABEL `optimal`, cross-platform identical; MOSEK
  ## independently agrees at 0.0164958). Fresh Problem object -- reused
  ## Constraint objects leak duals across solver calls via `.cache`.
  prob2 <- Problem(Minimize(matrix_frac(c(0, 1, 1), B)),
                   list(w >= 0, sum_entries(w) == 1))
  val2 <- psolve(prob2, solver = "CLARABEL")
  expect_equal(status(prob2), "optimal")
  expect_equal(as.numeric(val2), 0.0164945, tolerance = 1e-3)
})
