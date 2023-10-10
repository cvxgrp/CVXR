context("test-g05-conic_solvers-all")
TOL <- 1e-5

a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test installed solvers", {
  # Test the list of installed solvers.
  prob <- Problem(Minimize(norm1(x) + 1.0), list(x == 0))
  for(solver in names(SOLVER_MAP_CONIC)) {
    if(solver %in% INSTALLED_SOLVERS) {
      result <- solve(prob, solver = solver)
      expect_equal(result$value, 1.0, tolerance = TOL)
      expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = TOL)
    } else {
      expect_error(solve(prob, solver = solver), 
                   paste("The solver", solver, "is not installed"), fixed = TRUE)
    }
  }
  
  for(solver in names(SOLVER_MAP_QP)) {
    if(solver %in% INSTALLED_SOLVERS) {
      result <- solve(prob, solver = solver)
      expect_equal(result$getValue(x), matrix(c(0,0)), tolerance = TOL)
    } else {
      expect_error(solve(prob, solver = solver), 
                   paste("The solver", solver, "is not installed"), fixed = TRUE)
    }
  }
})

test_that("test mixed integer behavior", {
  x <- Variable(2, name = "x", integer = TRUE)
  objective <- Minimize(sum(x))
  prob <- Problem(objective, list(x >= 0))
  if(identical(INSTALLED_MI_SOLVERS, list("ECOS_BB"))) {
    expect_error(solve(prob), "You need a mixed-integer solver for this model")
  } else {
    result <- solve(prob)
    expect_equal(result$getValue(x), matrix(c(0,0)))
  }
})
