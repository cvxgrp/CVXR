context("test-g01-nonlinear_atoms")
TOL <- 1e-6

x <- Variable(2, name = "x")
y <- Variable(2, name = "y")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("Test log problem", {
  skip_on_cran()
  # Log in objective
  obj <- Maximize(sum(log(x)))
  constr <- list(x <= as.matrix(c(1, exp(1))))
  p <- Problem(obj, constr)
  result <- solve(p)
  expect_equal(result$value, 1, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(1, exp(1))), tolerance = TOL)

  # Log in constraint
  obj <- Minimize(sum(x))
  constr <- list(log(x) >= 0, x <= as.matrix(c(1,1)))
  p <- Problem(obj, constr)
  result <- solve(p)
  expect_equal(result$value, 2, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(1, 1)), tolerance = TOL)

  # Index into log
  obj <- Maximize(log(x)[2])
  constr <- list(x <= as.matrix(c(1, exp(1))))
  p <- Problem(obj, constr)
  result <- solve(p)
  expect_equal(result$value, 1, tolerance = TOL)

  # Scalar log
  obj <- Maximize(log(x[2]))
  constr <- list(x <= as.matrix(c(1, exp(1))))
  p <- Problem(obj, constr)
  result <- solve(p)
  expect_equal(result$value, 1, tolerance = TOL)
})

test_that("Test the entr function", {
  skip_on_cran()
  expect_equal(value(entr(0)), 0)
  expect_warning(expect_equal(value(entr(-1)), -Inf))
})

test_that("Test a problem with KL-divergence", {
  skip_on_cran()
  kK <- 50
  kSeed <- 10

  # Generate a random reference distribution
  # set.seed(kSeed)
  npSPriors <- matrix(stats::runif(kK), nrow = kK, ncol = 1)
  npSPriors <- npSPriors/sum(npSPriors)

  # Reference distribution
  p_refProb <- Parameter(kK, 1, nonneg = TRUE)

  # Distribution to be estimated
  v_prob <- Variable(kK, 1)
  objkl <- 0.0
  con <- 0.0
  for(k in 1:kK) {
    # objkl <- objkl + kl_div(v_prob[k,1], p_refProb[k,1])   # TODO: Parameters are unimplemented
    objkl <- objkl + kl_div(v_prob[k,1], npSPriors[k,1])
    con <- con + v_prob[k,1]
  }

  constrs <- list(con == 1)
  klprob <- Problem(Minimize(objkl), constrs)
  value(p_refProb) <- npSPriors

  result <- solve(klprob, solver = "SCS", verbose = TRUE)
  expect_equal(result$getValue(v_prob), npSPriors, tolerance = 1e-3)
  result <- solve(klprob, solver = "ECOS", verbose = TRUE)
  expect_equal(result$getValue(v_prob), npSPriors, tolerance = 1e-3)
})

test_that("Test a problem with entr", {
  skip_on_cran()
  for(n in c(5, 10, 25)) {
    print(n)
    x <- Variable(n)
    obj <- Maximize(sum(entr(x)))
    p <- Problem(obj, list(sum(x) == 1))
    result <- solve(p, solver = "ECOS", verbose = TRUE)
    expect_equal(result$getValue(x), matrix(rep(1.0/n, n)), tolerance = TOL)
    result <- solve(p, solver = "SCS", verbose = TRUE)
    expect_equal(result$getValue(x), matrix(rep(1.0/n, n)), tolerance = 1e-3)
  }
})

test_that("Test a problem with exp", {
  skip_on_cran()
  for(n in c(5, 10, 25)) {
    print(n)
    x <- Variable(n)
    obj <- Minimize(sum(exp(x)))
    p <- Problem(obj, list(sum(x) == 1))
    result <- solve(p, solver = "ECOS", verbose = TRUE)
    expect_equal(result$getValue(x), matrix(rep(1.0/n, n)), tolerance = TOL)
    result <- solve(p, solver = "SCS", verbose = TRUE)
    expect_equal(result$getValue(x), matrix(rep(1.0/n, n)), tolerance = 1e-3)
  }
})

test_that("Test a problem with log", {
  skip_on_cran()
  for(n in c(5, 10, 25)) {
    print(n)
    x <- Variable(n)
    obj <- Maximize(sum(log(x)))
    p <- Problem(obj, list(sum(x) == 1))
    result <- solve(p, solver = "ECOS", verbose = TRUE)
    expect_equal(result$getValue(x), matrix(rep(1.0/n, n)), tolerance = TOL)
    result <- solve(p, solver = "SCS", verbose = TRUE)
    expect_equal(result$getValue(x), matrix(rep(1.0/n, n)), tolerance = 1e-2)
  }
})
