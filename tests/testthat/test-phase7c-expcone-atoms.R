## Phase 7c: ExpCone + AxisAtom Atoms Tests
## Tests: kl_div, rel_entr, log1p, xexp, logistic, log_sum_exp, total_variation

# ═══════════════════════════════════════════════════════════════════
# KL Divergence: kl_div(x, y) = x*log(x/y) - x + y
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("KlDiv class construction and shape", {
  x <- Variable(c(3L, 1L))
  y <- Variable(c(3L, 1L))
  e <- kl_div(x, y)
  expect_s3_class(e, "CVXR::KlDiv")
  expect_equal(e@shape, c(3L, 1L))
  expect_length(e@args, 2L)
})

## @cvxpy NONE
test_that("KlDiv DCP properties match CVXPY", {
  x <- Variable(2, nonneg = TRUE)
  y <- Variable(2, nonneg = TRUE)
  e <- kl_div(x, y)

  ## Sign: always nonneg
  s <- sign_from_args(e)
  expect_true(s$is_nonneg)
  expect_false(s$is_nonpos)

  ## Curvature: convex
  expect_true(is_atom_convex(e))
  expect_false(is_atom_concave(e))

  ## Monotonicity: not monotone in either
  expect_false(is_incr(e, 1L))
  expect_false(is_incr(e, 2L))
  expect_false(is_decr(e, 1L))
  expect_false(is_decr(e, 2L))
})

## @cvxpy NONE
test_that("KlDiv domain constraints", {
  x <- Variable(2)
  y <- Variable(2)
  e <- kl_div(x, y)
  dom <- atom_domain(e)
  expect_length(dom, 2L)  # x >= 0, y >= 0
})

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_kl_div
test_that("KlDiv numeric values match scipy.special.kl_div", {
  e <- KlDiv(Constant(1), Constant(1))
  expect_equal(as.numeric(value(e)), 0.0, tolerance = 1e-10)

  e <- KlDiv(Constant(2), Constant(1))
  expect_equal(as.numeric(value(e)), 0.3862943611198906, tolerance = 1e-10)

  e <- KlDiv(Constant(0), Constant(1))
  expect_equal(as.numeric(value(e)), 1.0, tolerance = 1e-10)

  e <- KlDiv(Constant(1), Constant(2))
  expect_equal(as.numeric(value(e)), 0.3068528194400546, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("KlDiv broadcasting works", {
  ## Scalar + vector
  e <- kl_div(Constant(1), Constant(c(1, 2)))
  expect_equal(e@shape, c(2L, 1L))
})

# ═══════════════════════════════════════════════════════════════════
# Relative Entropy: rel_entr(x, y) = x*log(x/y)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("RelEntr class construction and shape", {
  x <- Variable(c(3L, 1L))
  y <- Variable(c(3L, 1L))
  e <- rel_entr(x, y)
  expect_s3_class(e, "CVXR::RelEntr")
  expect_equal(e@shape, c(3L, 1L))
})

## @cvxpy NONE
test_that("RelEntr DCP properties match CVXPY", {
  x <- Variable(2, nonneg = TRUE)
  y <- Variable(2, nonneg = TRUE)
  e <- rel_entr(x, y)

  ## Sign: unknown
  s <- sign_from_args(e)
  expect_false(s$is_nonneg)
  expect_false(s$is_nonpos)

  ## Curvature: convex
  expect_true(is_atom_convex(e))
  expect_false(is_atom_concave(e))

  ## Monotonicity: not incr in either; decr in arg 2 only
  expect_false(is_incr(e, 1L))
  expect_false(is_incr(e, 2L))
  expect_false(is_decr(e, 1L))
  expect_true(is_decr(e, 2L))
})

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_rel_entr
test_that("RelEntr numeric values match scipy.special.rel_entr", {
  expect_equal(as.numeric(value(RelEntr(Constant(1), Constant(1)))), 0.0, tolerance = 1e-10)
  expect_equal(as.numeric(value(RelEntr(Constant(2), Constant(1)))), 1.3862943611198906, tolerance = 1e-10)
  expect_equal(as.numeric(value(RelEntr(Constant(1), Constant(2)))), -0.6931471805599453, tolerance = 1e-10)
  expect_equal(as.numeric(value(RelEntr(Constant(0), Constant(1)))), 0.0, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════
# Log1p: log(1 + x)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Log1p class construction — extends Log", {
  x <- Variable(c(2L, 3L))
  e <- log1p_atom(x)
  expect_s3_class(e, "CVXR::Log1p")
  expect_s3_class(e, "CVXR::Log")  ## inherits from Log
  expect_equal(e@shape, c(2L, 3L))
})

## @cvxpy test_atoms.py::TestAtoms::test_log1p
test_that("Log1p DCP properties match CVXPY", {
  x <- Variable(2, nonneg = TRUE)
  e <- log1p_atom(x)

  ## Sign: tracks argument (nonneg arg → nonneg)
  s <- sign_from_args(e)
  expect_true(s$is_nonneg)
  expect_false(s$is_nonpos)

  ## Free variable → unknown sign
  z <- Variable(2)
  e2 <- log1p_atom(z)
  s2 <- sign_from_args(e2)
  expect_false(s2$is_nonneg)
  expect_false(s2$is_nonpos)

  ## Curvature: concave (inherited from Log)
  expect_false(is_atom_convex(e))
  expect_true(is_atom_concave(e))

  ## Monotonicity: increasing (inherited from Log)
  expect_true(is_incr(e, 1L))
  expect_false(is_decr(e, 1L))
})

## @cvxpy NONE
test_that("Log1p domain is x >= -1", {
  x <- Variable(2)
  e <- log1p_atom(x)
  dom <- atom_domain(e)
  expect_length(dom, 1L)
  ## The constraint should be x >= -1 (i.e., x + 1 >= 0)
})

## @cvxpy NONE
test_that("Log1p numeric values", {
  expect_equal(as.numeric(value(Log1p(Constant(0)))), 0.0, tolerance = 1e-10)
  expect_equal(as.numeric(value(Log1p(Constant(1)))), log(2), tolerance = 1e-10)
  expect_equal(as.numeric(value(Log1p(Constant(-0.5)))), log(0.5), tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════
# Xexp: x * exp(x)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Xexp class construction", {
  x <- Variable(c(2L, 1L))
  e <- xexp(x)
  expect_s3_class(e, "CVXR::Xexp")
  expect_equal(e@shape, c(2L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_xexp
test_that("Xexp DCP properties match CVXPY", {
  x <- Variable(2, nonneg = TRUE)
  e <- xexp(x)

  ## Sign: tracks arg (nonneg arg → nonneg)
  s <- sign_from_args(e)
  expect_true(s$is_nonneg)
  expect_false(s$is_nonpos)

  ## Curvature: convex when arg is nonneg
  expect_true(is_atom_convex(e))
  expect_false(is_atom_concave(e))

  ## Free variable → NOT convex (conditionally convex)
  z <- Variable(2)
  e2 <- xexp(z)
  expect_false(is_atom_convex(e2))

  ## Monotonicity: always increasing
  expect_true(is_incr(e, 1L))
  expect_false(is_decr(e, 1L))
})

## @cvxpy NONE
test_that("Xexp domain is x >= 0", {
  x <- Variable(2)
  e <- xexp(x)
  dom <- atom_domain(e)
  expect_length(dom, 1L)
})

## @cvxpy NONE
test_that("Xexp numeric values", {
  expect_equal(as.numeric(value(Xexp(Constant(0)))), 0.0, tolerance = 1e-10)
  expect_equal(as.numeric(value(Xexp(Constant(1)))), exp(1), tolerance = 1e-10)
  expect_equal(as.numeric(value(Xexp(Constant(2)))), 2 * exp(2), tolerance = 1e-10)
})

## @cvxpy NONE
test_that("Xexp log-log curvature", {
  x <- Variable(2, nonneg = TRUE)
  e <- xexp(x)
  expect_true(is_atom_log_log_convex(e))
  expect_false(is_atom_log_log_concave(e))
})

# ═══════════════════════════════════════════════════════════════════
# Logistic: log(1 + exp(x))
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Logistic class construction", {
  x <- Variable(c(3L, 2L))
  e <- logistic(x)
  expect_s3_class(e, "CVXR::Logistic")
  expect_equal(e@shape, c(3L, 2L))
})

## @cvxpy NONE
test_that("Logistic DCP properties match CVXPY", {
  z <- Variable(2)
  e <- logistic(z)

  ## Sign: always nonneg
  s <- sign_from_args(e)
  expect_true(s$is_nonneg)
  expect_false(s$is_nonpos)

  ## Curvature: convex
  expect_true(is_atom_convex(e))
  expect_false(is_atom_concave(e))

  ## Monotonicity: increasing
  expect_true(is_incr(e, 1L))
  expect_false(is_decr(e, 1L))
})

## @cvxpy NONE
test_that("Logistic numeric values (numerically stable)", {
  ## Using logaddexp(0, x) reference values
  expect_equal(as.numeric(value(Logistic(Constant(0)))), log(2), tolerance = 1e-10)
  expect_equal(as.numeric(value(Logistic(Constant(1)))), 1.3132616875182228, tolerance = 1e-10)
  expect_equal(as.numeric(value(Logistic(Constant(-1)))), 0.31326168751822286, tolerance = 1e-10)

  ## Numerical stability for large values
  expect_equal(as.numeric(value(Logistic(Constant(100)))), 100, tolerance = 1e-10)
  expect_equal(as.numeric(value(Logistic(Constant(-100)))), 0, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════
# LogSumExp: log(sum(exp(x)))
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("LogSumExp class construction — extends AxisAtom", {
  x <- Variable(c(3L, 1L))
  e <- log_sum_exp(x)
  expect_s3_class(e, "CVXR::LogSumExp")
  expect_s3_class(e, "CVXR::AxisAtom")
  expect_equal(e@shape, c(1L, 1L))  ## axis=NULL → scalar
})

## @cvxpy NONE
test_that("LogSumExp shapes with axis", {
  x <- Variable(c(3L, 4L))

  ## axis=NULL → (1, 1)
  expect_equal(log_sum_exp(x)@shape, c(1L, 1L))

  ## axis=2 → (1, 4) reduced column-wise
  expect_equal(log_sum_exp(x, axis = 2L)@shape, c(1L, 4L))

  ## axis=1 → (3, 1) reduced cols
  expect_equal(log_sum_exp(x, axis = 1L)@shape, c(3L, 1L))

  ## keepdims
  expect_equal(log_sum_exp(x, axis = 2L, keepdims = TRUE)@shape, c(1L, 4L))
  expect_equal(log_sum_exp(x, axis = 1L, keepdims = TRUE)@shape, c(3L, 1L))
})

## @cvxpy test_atoms.py::TestAtoms::test_log_sum_exp
test_that("LogSumExp DCP properties match CVXPY", {
  x <- Variable(3, nonneg = TRUE)
  e <- log_sum_exp(x)

  ## Sign: nonneg when arg nonneg
  s <- sign_from_args(e)
  expect_true(s$is_nonneg)
  expect_false(s$is_nonpos)

  ## Free variable → unknown sign
  z <- Variable(3)
  e2 <- log_sum_exp(z)
  s2 <- sign_from_args(e2)
  expect_false(s2$is_nonneg)

  ## Curvature: convex
  expect_true(is_atom_convex(e))
  expect_false(is_atom_concave(e))

  ## Monotonicity: increasing
  expect_true(is_incr(e, 1L))
  expect_false(is_decr(e, 1L))
})

## @cvxpy NONE
test_that("LogSumExp numeric values match scipy.special.logsumexp", {
  ## axis=NULL
  e <- LogSumExp(Constant(c(1, 2, 3)))
  expect_equal(as.numeric(value(e)), 3.40760596444438, tolerance = 1e-10)

  ## axis=0
  m <- matrix(c(1, 3, 2, 4), nrow = 2, ncol = 2)  # col-major: [[1,2],[3,4]]
  e0 <- LogSumExp(Constant(m), axis = 2L)
  expected <- c(3.12692801, 4.12692801)
  expect_equal(as.numeric(value(e0)), expected, tolerance = 1e-5)

  ## axis=1
  e1 <- LogSumExp(Constant(m), axis = 1L)
  expected <- c(2.31326169, 4.31326169)
  expect_equal(as.numeric(value(e1)), expected, tolerance = 1e-5)
})

# ═══════════════════════════════════════════════════════════════════
# Total Variation
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("total_variation 1D vector", {
  ## tv of a constant vector: |diff(x)|_1
  v <- c(1, 3, 2, 5)
  expect_equal(as.numeric(value(total_variation(Constant(v)))),
               sum(abs(diff(v))), tolerance = 1e-10)
})

## @cvxpy NONE
test_that("total_variation errors on scalar", {
  expect_error(total_variation(Constant(5)), "scalar")
})

## @cvxpy NONE
test_that("total_variation 2D matrix returns Expression", {
  x <- Variable(c(3L, 4L))
  e <- total_variation(x)
  expect_true(S7_inherits(e, Expression))
})

# ═══════════════════════════════════════════════════════════════════
# CANON_METHODS Registration
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("DCP canonicalization S7 methods registered for exp/power cone atoms", {
  x <- Variable()
  y <- Variable(nonneg = TRUE)
  expect_true(has_dcp_canon(KlDiv(y, y)))
  expect_true(has_dcp_canon(RelEntr(y, y)))
  expect_true(has_dcp_canon(Log1p(y)))
  expect_true(has_dcp_canon(Xexp(y)))
  expect_true(has_dcp_canon(Logistic(x)))
  expect_true(has_dcp_canon(LogSumExp(Variable(3))))
})

# ═══════════════════════════════════════════════════════════════════
# Canonicalization Constraint Counts
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("kl_div_canon produces 1 ExpCone constraint", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(2L, 1L))
  e <- kl_div(x, y)
  result <- kl_div_canon(e, e@args)
  expect_length(result[[2L]], 1L)
  expect_s3_class(result[[2L]][[1L]], "CVXR::ExpCone")
})

## @cvxpy NONE
test_that("rel_entr_canon produces 1 ExpCone constraint", {
  x <- Variable(c(2L, 1L))
  y <- Variable(c(2L, 1L))
  e <- rel_entr(x, y)
  result <- rel_entr_canon(e, e@args)
  expect_length(result[[2L]], 1L)
  expect_s3_class(result[[2L]][[1L]], "CVXR::ExpCone")
})

## @cvxpy NONE
test_that("logistic_canon produces correct constraints (2 ExpCone + 1 ineq)", {
  x <- Variable(c(2L, 1L))
  e <- logistic(x)
  result <- logistic_canon(e, e@args)
  obj <- result[[1L]]
  constrs <- result[[2L]]
  ## 2 ExpCone from exp_canon + 1 inequality (t1 + t2 <= 1)
  expect_equal(length(constrs), 3L)
  expect_s3_class(constrs[[1L]], "CVXR::ExpCone")
  expect_s3_class(constrs[[2L]], "CVXR::ExpCone")
})

## @cvxpy NONE
test_that("log1p_canon delegates to log_canon", {
  x <- Variable(c(2L, 1L))
  e <- log1p_atom(x)
  result <- log1p_canon(e, e@args)
  ## Should produce 1 ExpCone (same as log_canon)
  expect_length(result[[2L]], 1L)
  expect_s3_class(result[[2L]][[1L]], "CVXR::ExpCone")
})

# ═══════════════════════════════════════════════════════════════════
# DCP Composition Tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("kl_div is DCP valid", {
  x <- Variable(2, nonneg = TRUE)
  y <- Variable(2, nonneg = TRUE)
  e <- sum_entries(kl_div(x, y))
  expect_true(is_dcp(e))
})

## @cvxpy NONE
test_that("rel_entr is DCP valid in minimization", {
  x <- Variable(2, nonneg = TRUE)
  e <- sum_entries(rel_entr(x, Constant(c(1, 2))))
  expect_true(is_dcp(e))
})

## @cvxpy NONE
test_that("logistic is DCP valid", {
  x <- Variable(2)
  e <- sum_entries(logistic(x))
  expect_true(is_dcp(e))
})

## @cvxpy NONE
test_that("log_sum_exp is DCP valid", {
  x <- Variable(3)
  e <- log_sum_exp(x)
  expect_true(is_dcp(e))
})

## @cvxpy NONE
test_that("log1p is DCP valid in maximization", {
  x <- Variable(2, nonneg = TRUE)
  prob <- Problem(Maximize(sum_entries(log1p_atom(x))), list(x <= 5))
  expect_true(is_dcp(prob))
})

## @cvxpy NONE
test_that("xexp is DCP valid when arg is nonneg", {
  x <- Variable(2, nonneg = TRUE)
  e <- sum_entries(xexp(x))
  expect_true(is_dcp(e))
})

# ═══════════════════════════════════════════════════════════════════
# Solve Parity vs CVXPY
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("kl_div solve matches CVXPY", {
  skip_if_not_installed("clarabel")
  x <- Variable(2, nonneg = TRUE)
  prob <- Problem(Minimize(sum_entries(kl_div(x, Constant(c(1, 2))))),
                  list(x >= 0.5))
  result <- psolve(prob, solver = "CLARABEL", verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## CVXPY: value ≈ 0 (x converges to [1, 2])
  expect_equal(result, 0, tolerance = 1e-4)
  xval <- as.numeric(value(x))
  expect_equal(xval, c(1, 2), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("rel_entr solve matches CVXPY", {
  skip_if_not_installed("clarabel")
  x <- Variable(2, nonneg = TRUE)
  prob <- Problem(Minimize(sum_entries(rel_entr(x, Constant(c(1, 2))))),
                  list(x >= 0.5))
  result <- psolve(prob, solver = "CLARABEL", verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## CVXPY: value ≈ -1.0823
  expect_equal(result, -1.0823324721844125, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("logistic solve matches CVXPY", {
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(logistic(x))),
                  list(x >= -1, x <= 1))
  result <- psolve(prob, solver = "CLARABEL", verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## CVXPY: value ≈ 0.6265, x ≈ [-1, -1]
  expect_equal(result, 0.6265233760210392, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(-1, -1), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("log_sum_exp solve matches CVXPY", {
  skip_if_not_installed("clarabel")
  x <- Variable(3)
  prob <- Problem(Minimize(log_sum_exp(x)),
                  list(x >= 0, sum_entries(x) == 3))
  result <- psolve(prob, solver = "CLARABEL", verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## CVXPY: value ≈ 2.0986, x ≈ [1, 1, 1]
  expect_equal(result, 2.0986122880312275, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(1, 1, 1), tolerance = 1e-2)
})

## @cvxpy NONE
test_that("log1p solve matches CVXPY", {
  skip_if_not_installed("clarabel")
  x <- Variable(2, nonneg = TRUE)
  prob <- Problem(Maximize(sum_entries(log1p_atom(x))),
                  list(x <= 2))
  result <- psolve(prob, solver = "CLARABEL", verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## CVXPY: value ≈ 2.1972 = 2*log(3)
  expect_equal(result, 2 * log(3), tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(2, 2), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("xexp solve matches CVXPY", {
  skip_if_not_installed("clarabel")
  x <- Variable(2, nonneg = TRUE)
  prob <- Problem(Minimize(sum_entries(xexp(x))),
                  list(x >= 0.5, x <= 2))
  result <- psolve(prob, solver = "CLARABEL", verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## CVXPY: value ≈ 1.6487, x ≈ [0.5, 0.5]
  expect_equal(result, 1.6487212704396041, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(0.5, 0.5), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("total_variation solve matches CVXPY", {
  skip_if_not_installed("clarabel")
  x <- Variable(5)
  prob <- Problem(Minimize(total_variation(x)),
                  list(x[1] == 0, x[5] == 4))
  result <- psolve(prob, solver = "CLARABEL", verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## CVXPY: value ≈ 4 (evenly spaced)
  expect_equal(result, 4.0, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════
# log_sum_exp with axis solve parity
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("log_sum_exp axis=0 solve", {
  skip_if_not_installed("clarabel")
  x <- Variable(c(3L, 2L))
  ## Minimize sum of column-wise log_sum_exp
  prob <- Problem(Minimize(sum_entries(log_sum_exp(x, axis = 2L))),
                  list(x >= 0, sum_entries(x) == 6))
  result <- psolve(prob, solver = "CLARABEL", verbose = FALSE)
  expect_equal(status(prob), "optimal")
  ## Optimal: uniform x = 1
  expect_equal(as.numeric(value(x)), rep(1, 6), tolerance = 0.1)
})
