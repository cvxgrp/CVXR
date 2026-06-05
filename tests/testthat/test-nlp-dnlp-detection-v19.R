## CVXPY 1.9 DNLP detection (nlp_tests/test_dnlp.py).
##
## Predicate tests for is_linearizable_convex / is_linearizable_concave /
## is_dnlp on smooth and piecewise-linear compositions, plus problems that are
## not DNLP and must raise a DNLP error under psolve(nlp = TRUE). All predicate
## values cross-checked against CVXPY 1.9.0. CVXPY's `minimum`/`maximum` map to
## CVXR's `min_elemwise`/`max_elemwise`; a constraint's residual expression is
## `constraint@.expr`.

## @cvxpy nlp_tests/test_dnlp.py::TestDNLP::test_abs_linearizable_convex
test_that("DNLP detection: abs is linearizable-convex and DNLP", {
  x <- Variable()
  expect_true(is_linearizable_convex(abs(x)))
  expect_true(is_dnlp(Problem(Minimize(abs(x)))))
})

## @cvxpy nlp_tests/test_dnlp.py::TestDNLP::test_sqrt_linearizable_concave
test_that("DNLP detection: sqrt is linearizable-concave and DNLP", {
  x <- Variable()
  expect_true(is_linearizable_concave(sqrt(x)))
  expect_true(is_dnlp(Problem(Maximize(sqrt(x)))))
})

## @cvxpy nlp_tests/test_dnlp.py::TestDNLP::test_simple_neg_expr
test_that("DNLP detection: abs(x) - sqrt(y) <= 5 is a DNLP constraint", {
  x <- Variable(); y <- Variable()
  con <- (abs(x) - sqrt(y) <= 5)
  expect_true(is_dnlp(con))
  expect_true(is_linearizable_convex(con@.expr))
})

## @cvxpy nlp_tests/test_dnlp.py::TestDNLP::test_non_dnlp
test_that("DNLP detection: abs(x) >= 5 is concave but not a DNLP constraint", {
  x <- Variable()
  con <- (abs(x) >= 5)
  expect_true(is_linearizable_concave(con@.expr))
  expect_false(is_dnlp(con))
})

## @cvxpy nlp_tests/test_dnlp.py::TestDNLP::test_simple_composition
test_that("DNLP detection: simple smooth-over-PWL compositions are DNLP", {
  x <- Variable()
  expect_true(is_dnlp(Minimize(log(abs(x)))))
  expect_true(is_dnlp(Minimize(exp(norm1(x)))))
  expect_true(is_dnlp(sqrt(abs(x))))
})

## @cvxpy nlp_tests/test_dnlp.py::TestDNLP::test_complicated_composition
test_that("DNLP detection: minimum(sqrt(exp(x)), -abs(y)) is concave, not DNLP to minimize", {
  x <- Variable(); y <- Variable()
  expr <- min_elemwise(sqrt(exp(x)), -abs(y))
  expect_true(is_linearizable_concave(expr))
  expect_false(is_dnlp(Problem(Minimize(expr))))
})

## @cvxpy nlp_tests/test_dnlp.py::TestNonDNLP::test_max
test_that("DNLP detection: maximize maximum(x, y) errors under nlp = TRUE", {
  x <- Variable(1); y <- Variable(1)
  prob <- Problem(Maximize(max_elemwise(x, y)), list(x - 14 == 0, y - 6 == 0))
  expect_error(psolve(prob, nlp = TRUE), "not DNLP")
})

## @cvxpy nlp_tests/test_dnlp.py::TestNonDNLP::test_min
test_that("DNLP detection: minimize minimum(x, y) errors under nlp = TRUE", {
  x <- Variable(1); y <- Variable(1)
  prob <- Problem(Minimize(min_elemwise(x, y)), list(x - 14 == 0, y - 6 == 0))
  expect_error(psolve(prob, nlp = TRUE), "not DNLP")
})

## @cvxpy nlp_tests/test_dnlp.py::TestNonDNLP::test_max_2
test_that("DNLP detection: maximize sum(maximum(x, y)) errors under nlp = TRUE", {
  x <- Variable(3); y <- Variable(3)
  prob <- Problem(Maximize(sum_entries(max_elemwise(x, y))), list(x <= 14, y <= 14))
  expect_error(psolve(prob, nlp = TRUE), "not DNLP")
})
