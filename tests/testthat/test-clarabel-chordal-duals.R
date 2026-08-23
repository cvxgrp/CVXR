## CLARABEL chordal decomposition must not change PSD dual variables.
##
## R clarabel 0.11.3 turned `chordal_decomposition_enable` on by default for
## semidefinite programs, and warned that "the dual variables reported for
## semidefinite constraints are not in general the same as before". For CVXR
## they must stay the same to solver tolerance, because 0.11.3 also enables
## `chordal_decomposition_complete_dual`, which reconstructs the full dual.
##
## Both settings are passed EXPLICITLY here rather than relying on the package
## default, so the test states the invariant regardless of which clarabel
## version is installed.

## @cvxpy NONE -- R clarabel packaging behavior, no CVXPY counterpart.
test_that("chordal decomposition leaves PSD duals unchanged", {
  skip_if_not_installed("clarabel")

  ## A PSD block with STRUCTURAL zeros off the band: M = tridiag(1, d, 1).
  ## A free symmetric `X >> 0` would NOT do -- every entry is its own
  ## variable, so there is no sparsity for the decomposition to exploit and
  ## the test would silently compare two identical code paths.
  n <- 20L
  d <- Variable(n)
  terms <- vector("list", n + 1L)
  for (i in seq_len(n)) {
    B <- matrix(0, n, n); B[i, i] <- 1
    terms[[i]] <- Constant(B) * d[i]
  }
  Boff <- matrix(0, n, n)
  for (i in seq_len(n - 1L)) { Boff[i, i + 1L] <- 1; Boff[i + 1L, i] <- 1 }
  terms[[n + 1L]] <- Constant(Boff)
  M <- Reduce(`+`, terms)

  solve_with <- function(chordal) {
    con <- M %>>% 0
    prob <- Problem(Minimize(sum(d)), list(con))
    val <- psolve(prob, solver = "CLARABEL",
                  chordal_decomposition_enable = chordal)
    list(val = val, dual = as.matrix(dual_value(con)), status = status(prob))
  }

  on  <- solve_with(TRUE)
  off <- solve_with(FALSE)

  expect_equal(on$status, "optimal")
  expect_equal(off$status, "optimal")
  expect_equal(on$val, off$val, tolerance = 1e-5)

  ## Same shape, and agreeing entrywise to solver tolerance.
  expect_equal(dim(on$dual), dim(off$dual))
  expect_equal(dim(on$dual), c(n, n))
  expect_lt(max(abs(on$dual - off$dual)), 1e-4)

  ## And the reconstructed dual must still be a valid one: PSD.
  for (S in list(on$dual, off$dual)) {
    Ssym <- (S + t(S)) / 2
    expect_gt(min(.eigvalsh(Ssym, only_values = TRUE)$values), -1e-6)
  }
})
