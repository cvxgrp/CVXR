atoms <- list(
  list(
    list(
      list(Abs, c(2,2), rbind(c(-5,2), c(-3,1)), Constant(rbind(c(5,2), c(3,1)))),
      list(Diag, c(2,1), rbind(c(-5,2), c(-3,1)), Constant(c(-5,1))),
      list(Diag, c(2,2), c(-5,1), Constant(rbind(c(-5,0), c(0,1)))),
      list(Exp, c(2,2), rbind(c(1,0), c(2,-1)), Constant(rbind(c(exp(1),1), c(exp(2), exp(-1))))),
      list(Huber, c(2,2), rbind(c(0.5,-1.5), c(4,0)), Constant(rbind(c(0.25,2), c(7,0)))),
        ),
    Minimize),
  list(
    list(
      list(Entr, c(2,2), rbind(c(1,exp(1)), c(exp(2), exp(-1))), Constant(rbind(c(0,-exp(1)), c(-2*exp(2),exp(-1)))))
    ),
    Maximize
  )
)

check_solver <- function(prob, solver_name) {
  canon <- canonicalize(prob)
  objective <- canon[[1]]
  constraints <- canon[[2]]
  validate_solver(solver_name, constraints)
}

run_atom <- function(atom, problem, obj_val, solver, verbose = FALSE) {
  expect_true(is_dcp(problem))
  if(verbose) {
    print(problem@objective)
    print(problem@constraints)
    print(paste("solver", solver))
  }
}
