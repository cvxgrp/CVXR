## CVXPY 1.9 test_scalarize.py parity -- the `scalarize` multi-objective transform.
##
## CVXPY's `cvxpy.transforms.scalarize` submodule maps to CVXR's exported
## `scalarize` list: scalarize$weighted_sum / $targets_and_priorities / $max /
## $log_sum_exp. Maximize objectives are built as Maximize(-obj@args[[1L]]),
## mirroring CVXPY's `cp.Maximize(-obj.args[0])`.

## Shared setUp: x scalar; obj1 = Minimize(x^2), obj2 = Minimize((x-1)^2).
.scz_setup <- function() {
  x <- Variable()
  list(x = x,
       objectives = list(Minimize(square(x)), Minimize(square(x - 1))))
}

## @cvxpy test_scalarize.py::ScalarizeTest::test_weighted_sum
test_that("weighted_sum combines objectives by weight", {
  s <- .scz_setup(); x <- s$x; objs <- s$objectives

  prob <- Problem(scalarize$weighted_sum(objs, c(1, 1)))
  psolve(prob)
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-3)

  prob <- Problem(scalarize$weighted_sum(objs, c(1, 0)))
  psolve(prob)
  expect_equal(as.numeric(value(x)), 0, tolerance = 1e-3)

  prob <- Problem(scalarize$weighted_sum(objs, c(0, 1)))
  psolve(prob)
  expect_equal(as.numeric(value(x)), 1, tolerance = 1e-3)
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_targets_and_priorities
test_that("targets_and_priorities penalizes within target/limit range", {
  s <- .scz_setup(); x <- s$x; objs <- s$objectives

  prob <- Problem(scalarize$targets_and_priorities(objs, c(1, 1), c(1, 1)))
  psolve(prob)
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-3)

  prob <- Problem(scalarize$targets_and_priorities(objs, c(1, 1), c(1, 0)))
  psolve(prob)
  expect_equal(as.numeric(value(x)), 1, tolerance = 1e-3)

  scz <- scalarize$targets_and_priorities(objs, c(1, 1e-4), c(0, 0),
                                          limits = c(1, 0.25), off_target = 1e-5)
  prob <- Problem(scz)
  psolve(prob)
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-3)

  ## Maximize variants
  max_objs <- lapply(objs, function(o) Maximize(-o@args[[1L]]))
  scz <- scalarize$targets_and_priorities(max_objs, c(1, 1), c(-1, 0), off_target = 1e-5)
  expect_true(is_concave(scz@args[[1L]]))
  prob <- Problem(scz)
  psolve(prob)
  expect_equal(as.numeric(value(x)), 1, tolerance = 1e-3)

  scz <- scalarize$targets_and_priorities(max_objs, c(1, 1e-4), c(0, 0),
                                          limits = c(-1, -0.25), off_target = 1e-5)
  expect_true(is_concave(scz@args[[1L]]))
  prob <- Problem(scz)
  psolve(prob)
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-3)
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_negative_priority_regression
test_that("targets_and_priorities uses flipped target/limit for negative priority", {
  s <- .scz_setup(); x <- s$x; objs <- s$objectives
  obj2 <- Maximize(-objs[[2L]]@args[[1L]])
  scz <- scalarize$targets_and_priorities(
    list(objs[[1L]], obj2), c(1, -1), c(1, -1),
    limits = c(0.5, -0.5), off_target = 1e-2)
  prob <- Problem(scz)
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-3)
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_mixed_convexity
test_that("targets_and_priorities resolves mixed convexity by curvature", {
  s <- .scz_setup(); objs <- s$objectives
  obj1 <- objs[[1L]]
  obj2 <- Maximize(-objs[[2L]]@args[[1L]])
  mixed <- list(obj1, obj2)

  expect_error(
    scalarize$targets_and_priorities(mixed, c(1, 1), c(1, -1)),
    "neither convex nor concave")

  scz <- scalarize$targets_and_priorities(mixed, c(1, -1), c(1, -1), limits = c(1, -1))
  expect_true(is_convex(scz@args[[1L]]))

  scz <- scalarize$targets_and_priorities(mixed, c(-1, 1), c(1, -1), limits = c(1, -1))
  expect_true(is_concave(scz@args[[1L]]))
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_targets_and_priorities_exceptions
test_that("targets_and_priorities validates argument lengths", {
  s <- .scz_setup(); objs <- s$objectives

  expect_error(scalarize$targets_and_priorities(objs, c(1), c(1, 1)),
               "Number of objectives and priorities")
  expect_error(scalarize$targets_and_priorities(objs, c(1, 1), c(1)),
               "Number of objectives and targets")
  expect_error(scalarize$targets_and_priorities(objs, c(1, 1), c(1, 1), limits = c(1)),
               "Number of objectives and limits")
  expect_error(
    scalarize$targets_and_priorities(objs, c(1, 1), c(1, 1), limits = c(1, 1),
                                     off_target = -1),
    "off_target argument must be nonnegative")
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_maximize_affine_targets_and_priorities
test_that("Maximize(affine) is not silently flipped to Minimize", {
  x <- Variable()
  obj <- Maximize(x)
  scz <- scalarize$targets_and_priorities(list(obj), 1, 5, off_target = 0.01)
  expect_s7_class(scz, Maximize)
  prob <- Problem(scz, list(x <= 10, x >= 0))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_true(as.numeric(value(x)) > 5.0)
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_minimize_affine_negative_priority
test_that("Minimize(affine) with negative priority flips to Maximize", {
  x <- Variable()
  obj <- Minimize(x)
  scz <- scalarize$targets_and_priorities(list(obj), -1, 5, off_target = 0.01)
  expect_s7_class(scz, Maximize)
  prob <- Problem(scz, list(x <= 10, x >= 0))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_true(as.numeric(value(x)) < 5.0)
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_maximize_priority_leq_off_target
test_that("Maximize with priority <= off_target stays Maximize / errors", {
  x <- Variable()

  scz <- scalarize$targets_and_priorities(list(Maximize(x)), 1e-5, 3, off_target = 1e-5)
  expect_s7_class(scz, Maximize)
  prob <- Problem(scz, list(x >= 0, x <= 10))
  psolve(prob, solver = "CLARABEL")
  expect_true(as.numeric(value(x)) > 5.0)

  expect_error(
    scalarize$targets_and_priorities(list(Maximize(x)), 0, 3, off_target = 1e-5),
    "not concave")

  y <- Variable()
  scz <- scalarize$targets_and_priorities(
    list(Maximize(x), Maximize(y)), c(1e-5, 1e-5), c(3, 7), off_target = 1e-5)
  expect_s7_class(scz, Maximize)
  prob <- Problem(scz, list(x >= 0, x <= 10, y >= 0, y <= 10))
  psolve(prob, solver = "CLARABEL")
  expect_true(as.numeric(value(x)) > 5.0)
  expect_true(as.numeric(value(y)) > 5.0)
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_minimize_priority_leq_off_target
test_that("Minimize(affine) with priority <= off_target stays Minimize / errors", {
  x <- Variable()

  expect_error(
    scalarize$targets_and_priorities(list(Minimize(x)), 0, 3, off_target = 1e-5),
    "not convex")

  scz <- scalarize$targets_and_priorities(list(Minimize(x)), 1e-5, 3, off_target = 1e-5)
  expect_s7_class(scz, Minimize)
  prob <- Problem(scz, list(x >= 0, x <= 10))
  psolve(prob, solver = "CLARABEL")
  expect_true(as.numeric(value(x)) < 5.0)
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_max
test_that("max minimizes the largest weighted objective term", {
  s <- .scz_setup(); x <- s$x; objs <- s$objectives

  prob <- Problem(scalarize$max(objs, c(1, 2)))
  psolve(prob)
  expect_equal(as.numeric(value(x)), 0.5858, tolerance = 1e-3)

  prob <- Problem(scalarize$max(objs, c(2, 1)))
  psolve(prob)
  expect_equal(as.numeric(value(x)), 0.4142, tolerance = 1e-3)
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_log_sum_exp
test_that("log_sum_exp combines objectives smoothly", {
  s <- .scz_setup(); x <- s$x; objs <- s$objectives

  prob <- Problem(scalarize$log_sum_exp(objs, c(1, 2)))
  psolve(prob)
  expect_equal(as.numeric(value(x)), 0.6354, tolerance = 1e-3)

  prob <- Problem(scalarize$log_sum_exp(objs, c(2, 1)))
  psolve(prob)
  expect_equal(as.numeric(value(x)), 0.3646, tolerance = 1e-3)
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_weighted_sum_maximize_and_negative_weights
test_that("weighted_sum handles Maximize objectives and negative weights", {
  x <- Variable()
  objs <- list(Maximize(-square(x)), Maximize(-square(x - 1)))
  scz <- scalarize$weighted_sum(objs, c(1, 1))
  expect_s7_class(scz, Maximize)
  prob <- Problem(scz)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-3)

  ## Negative weight on a Minimize flips to Maximize.
  scz <- scalarize$weighted_sum(list(Minimize(square(x))), -1)
  expect_s7_class(scz, Maximize)

  ## Mixed types via negative weight -> DCP error on objective addition.
  s <- .scz_setup()
  expect_error(scalarize$weighted_sum(s$objectives, c(-1, 1)),
               "Cannot add")
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_targets_and_priorities_penalty_values_minimize
test_that("targets_and_priorities Minimize penalty matches the formula", {
  x <- Variable()
  scz <- scalarize$targets_and_priorities(list(Minimize(square(x))), 2.0, 1.0,
                                          off_target = 0.1)

  prob <- Problem(scz, list(x == 0)); psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 0.0, tolerance = 1e-4)

  prob <- Problem(scz, list(x == 0.5)); psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 0.025, tolerance = 1e-4)

  prob <- Problem(scz, list(x == 2)); psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 6.1, tolerance = 1e-3)
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_targets_and_priorities_penalty_values_maximize
test_that("targets_and_priorities Maximize penalty matches the formula", {
  x <- Variable()
  obj <- Maximize(-square(x))
  scz <- scalarize$targets_and_priorities(list(obj), 2.0, -1.0, off_target = 0.1)
  expect_s7_class(scz, Maximize)

  prob <- Problem(scz, list(x == 0)); psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), 0.0, tolerance = 1e-4)

  prob <- Problem(scz, list(x == 2)); psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(prob)), -6.1, tolerance = 1e-3)
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_maximize_negative_priority_with_limits
test_that("Maximize + negative priority + limits triple-flips to Minimize", {
  x <- Variable()
  obj <- Maximize(-square(x))
  scz <- scalarize$targets_and_priorities(list(obj), -1, -1, limits = -4,
                                          off_target = 0.01)
  expect_s7_class(scz, Minimize)
  prob <- Problem(scz, list(x >= -5, x <= 5))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 0.0, tolerance = 1e-2)
  expect_true(abs(as.numeric(value(x))) <= 2.0 + 1e-3)
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_limits_enforced_minimize_and_maximize
test_that("limits act as hard bounds for Minimize and Maximize", {
  x <- Variable()

  scz <- scalarize$targets_and_priorities(list(Minimize(x)), 1, 0, limits = 3,
                                          off_target = 0.01)
  prob <- Problem(scz, list(x >= -1, x <= 10)); psolve(prob, solver = "CLARABEL")
  expect_true(as.numeric(value(x)) <= 3.0 + 1e-3)
  expect_equal(as.numeric(value(x)), -1.0, tolerance = 1e-2)

  scz <- scalarize$targets_and_priorities(list(Maximize(x)), 1, 8, limits = 2,
                                          off_target = 0.01)
  prob <- Problem(scz, list(x >= 0, x <= 10)); psolve(prob, solver = "CLARABEL")
  expect_true(as.numeric(value(x)) >= 2.0 - 1e-3)
  expect_equal(as.numeric(value(x)), 10.0, tolerance = 1e-2)

  scz <- scalarize$targets_and_priorities(list(Minimize(x)), 1, 0, limits = -1,
                                          off_target = 0.01)
  prob <- Problem(scz, list(x >= 0, x <= 10)); psolve(prob, solver = "CLARABEL")
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_off_target_edge_cases
test_that("targets_and_priorities off_target edge cases", {
  x <- Variable()
  obj <- Minimize(square(x))

  scz <- scalarize$targets_and_priorities(list(obj), 1, 10, off_target = 0)
  prob <- Problem(scz, list(x >= -5, x <= 5)); psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 0.0, tolerance = 1e-2)

  scz <- scalarize$targets_and_priorities(list(obj), 0.5, 10, off_target = 0.5)
  prob <- Problem(scz); psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 0.0, tolerance = 1e-3)

  expect_error(
    scalarize$targets_and_priorities(list(obj), 0.1, 1, off_target = 1.0),
    "not convex")

  obj_aff <- Minimize(x)
  expect_error(
    scalarize$targets_and_priorities(list(obj_aff), 0.1, 1, off_target = 1.0),
    "not convex")
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_all_negative_priorities
test_that("all-negative priorities flip every objective", {
  x <- Variable()
  objs <- list(Minimize(square(x)), Minimize(square(x - 1)))
  scz <- scalarize$targets_and_priorities(objs, c(-1, -1), c(-1, -1), off_target = 1e-5)
  expect_s7_class(scz, Maximize)
  prob <- Problem(scz); psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-3)
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_log_sum_exp_gamma_convergence
test_that("log_sum_exp converges to weighted_sum (gamma->0) and max (gamma->Inf)", {
  s <- .scz_setup(); x <- s$x; objs <- s$objectives

  prob <- Problem(scalarize$weighted_sum(objs, c(1, 1)))
  psolve(prob, solver = "CLARABEL"); x_ws <- as.numeric(value(x))

  prob <- Problem(scalarize$max(objs, c(1, 1)))
  psolve(prob, solver = "CLARABEL"); x_max <- as.numeric(value(x))

  prob <- Problem(scalarize$log_sum_exp(objs, c(1, 1), gamma = 0.01))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), x_ws, tolerance = 1e-2)

  prob <- Problem(scalarize$log_sum_exp(objs, c(1, 1), gamma = 100))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), x_max, tolerance = 1e-2)
})

## @cvxpy test_scalarize.py::ScalarizeTest::test_max_and_log_sum_exp_maximize_not_dcp
test_that("max and log_sum_exp of Maximize objectives are not DCP", {
  x <- Variable()
  objs <- list(Maximize(-square(x)), Maximize(-square(x - 1)))

  expect_false(is_dcp(Problem(scalarize$max(objs, c(1, 1)))))
  expect_false(is_dcp(Problem(scalarize$log_sum_exp(objs, c(1, 1)))))
})
