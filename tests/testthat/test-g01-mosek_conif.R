context("test-g01-mosek_conif")

MOSEK_AVAILABLE  <- "MOSEK" %in% installed_solvers()
TOL <- 1e-6

if (MOSEK_AVAILABLE) library("Rmosek")

test_that("test MOSEK SOCP", {
  skip_on_cran()
  skip_if_not(MOSEK_AVAILABLE, "Skipping MOSEK test as it is not available.!")

  # Formulate the following SOCP with CVXR
  #    min 3 * x[1] + 2 * x[2] + x[3]
  #       s.t. p_norm(x,2) <= y[1]
  #            p_norm(x,2) <= y[2]
  #            x[1] + x[2] + 3*x[3] >= 1.0
  #            y <= 5
  # and solve with MOSEK and ECOS. Compare MOSEK and ECOS primal and dual solutions.

  x <- Variable(3)
  y <- Variable(2)
  constraints <- list(p_norm(x, 2) <= y[1],
                      p_norm(x, 2) <= y[2],
                      x[1] + x[2] + 3 * x[3] >= 1.0,
                      y <= 5)
  obj <- Minimize(3 * x[1] + 2 * x[2] + x[3])
  prob <- Problem(obj, constraints)
  result_mosek <- solve(prob, solver = "MOSEK")
  x_mosek <- result_mosek$getValue(x)
  duals_mosek <- lapply(constraints, function(c) { result_mosek$getDualValue(c) })

  result_ecos <- solve(prob, solver = "ECOS")
  x_ecos <- result_ecos$getValue(x)
  duals_ecos <- lapply(constraints, function(c) { result_ecos$getDualValue(c) })

  expect_equal(x_mosek, x_ecos, tolerance = TOL)
  expect_equal(length(duals_ecos), length(duals_mosek))
  for(i in seq_along(duals_mosek)) {
      expect_equal(duals_mosek[[i]], duals_ecos[[i]], tolerance = 1e-4)
  }

})

test_that("test MOSEK SDP", {
  skip_on_cran()
  skip_if_not(MOSEK_AVAILABLE, "Skipping MOSEK test as it is not available.!")
  # Solve "Example 8.3" from Convex Optimization by Boyd & Vandenberghe.
  # Verify (1) optimal objective values, (2) that the dual variable to the PSD constraint
  # belongs to the correct cone (i.e. the dual variable is itself PSD), and (3) that
  # complementary slackness holds with the PSD primal variable and its dual variable.

  # This is an example from Convex Optimization by B&V.
  # Example 8.3 (page 408 in the 19th printing).
  rho <- Variable(4, 4)
  constraints <- list(0.6 <= rho[1, 2], rho[1, 2] <= 0.9,
                      0.8 <= rho[1, 3], rho[1, 3] <= 0.9,
                      0.5 <= rho[2, 4], rho[2, 4] <= 0.7,
                      -0.8 <= rho[3, 4], rho[3, 4] <= -0.4,
                      rho[1, 1] == 1, rho[2, 2] == 1, rho[3, 3] == 1, rho[4, 4] == 1,
                      rho %>>% 0)

  # First, minimize rho[1, 4]
  obj <- Minimize(rho[1, 4])
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "MOSEK")
  expect_equal(round(result$value, 2), -0.39)
  Y <- result$getDualValue(constraints[[length(constraints)]])
  eigs <- eigen(Y, only.values = TRUE)
  expect_true(all(eigs$value >= 0))
  complementary_slackness <- sum(diag(result$getValue(rho) %*% Y))
  expect_equal(complementary_slackness, 0.0, tolereance = TOL)

  # Then, maximize rho[1, 4]
  obj <- Maximize(rho[1, 4])
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "MOSEK")
  expect_equal(round(result$value, 2), 0.23)
  Y <- result$getDualValue(constraints[[length(constraints)]])
  eigs <- eigen(Y, only.values = TRUE)
  expect_true(all(eigs$value >= 0))
  complementary_slackness <- sum(diag(result$getValue(rho) %*% Y))
  expect_equal(complementary_slackness, 0.0, tolerance = TOL)
})

test_that("test MOSEK exponential cone", {
  skip_on_cran()
  skip_if_not(MOSEK_AVAILABLE, "Skipping MOSEK test as it is not available.!")
  # Formulate the following exponential cone problem with cvxpy
  #   min   3 * x[1] + 2 * x[2] + x[3]
  #     s.t.  0.1 <= x[1] + x[2] + x[3] <= 1
  #           x >= 0
  #           x[1] >= x[2] * exp(x[3] / x[2])
  # and solve with MOSEK and ECOS. Ensure that MOSEK and ECOS have the same
  # primal and dual solutions.
  #
  # Note that the exponential cone constraint can be rewritten in terms of the
  # relative entropy cone. The correspondence is as follows:
  #     x[1] >= x[2] * exp(x[3] / x[2])
  #       iff
  #     x[2] * log(x[2] / x[1]) + x[3] <= 0.

  ## TODO: Only run if exponential cone is supported.
  ## Formulate and solve the problem with CVXR
  x <- Variable(3, 1)
  constraints <- list(sum(x) <= 1.0, sum(x) >= 0.1, x >= 0.01,
                      kl_div(x[2], x[1]) + x[2] - x[1] + x[3] <= 0)
  obj <- Minimize(3 * x[1] + 2 * x[2] + x[1])
  prob <- Problem(obj, constraints)
  result_mosek <- solve(prob, solver = "MOSEK")
  val_mosek = result_mosek$value
  x_mosek = result_mosek$getValue(x)
  duals_mosek = lapply(constraints, function(c) { result_mosek$getDualValue(c) })

  result_ecos <- solve(prob, solver = "ECOS")
  val_ecos = result_ecos$value
  x_ecos = result_ecos$getValue(x)
  duals_ecos = lapply(constraints, function(c) { result_ecos$getDualValue(c) })

  ## Verify results
  expect_equal(val_mosek, val_ecos, tolerance = TOL)
  expect_equal(x_mosek, x_ecos, tolerance = 1e-4)
  expect_equal(length(duals_ecos), length(duals_mosek))
  for(i in seq_along(duals_mosek))
      expect_equal(duals_mosek[[i]], duals_ecos[[i]], tolerance = 1e-4)
})

test_that("test MOSEK MI SOCP", {
  skip_on_cran()
  skip_if_not(MOSEK_AVAILABLE, "Skipping MOSEK test as it is not available.!")

  # Formulate the following mixed-integer SOCP with CVXR
  #   min 3 * x[1] + 2 * x[2] + x[3] + y[1] + 2 * y[2]
  #     s.t. p_norm(x,2) <= y[1]
  #          p_norm(x,2) <= y[2]
  #          x[1] + x[2] + 3*x[3] >= 0.1
  #          y <= 5, y integer.
  # and solve with MOSEK.

  x <- Variable(3)
  y <- Variable(2, integer = TRUE)
  constraints <- list(p_norm(x, 2) <= y[1],
                      p_norm(x, 2) <= y[2],
                      x[1] + x[2] + 3 * x[3] >= 0.1,
                      y <= 5)
  obj <- Minimize(3 * x[1] + 2 * x[2] + x[3] + y[1] + 2 * y[2])
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "MOSEK")
  expect_equal(result$getValue(y), as.matrix(c(1,1)), tolerance = TOL)
})

test_that("test MOSEK LP solution selection", {
  skip_on_cran()
  skip_if_not(MOSEK_AVAILABLE, "Skipping MOSEK test as it is not available.!")
  # Problem from
  # http://cvxopt.org/userguide/coneprog.html?highlight=solvers.lp#cvxopt.solvers.lp

  x <- Variable(2)
  objective <- Minimize(-4 * x[1] - 5 * x[2])
  constraints <- list(2 * x[1] + x[2] <= 3,
                      x[1] + 2 * x[2] <= 3,
                      x[1] >= 0, x[2] >= 0)
  prob <- Problem(objective, constraints)

  # Default solution (interior point)
  result <- solve(prob, solver = "MOSEK")
  expect_equal(result$value, -9, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(1,1)), tolerance = 1e-5)

  # NOTE: Rmosek does not provide the bfs option.
  # Basic feasible solution
  # result <- solve(prob, solver = "MOSEK", bfs = TRUE)
  # expect_equal(result$value, -9, tolerance = 1e-10)
  # expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = 1e-10)

})
