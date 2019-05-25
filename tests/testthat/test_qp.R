context("test-g01-qp")

a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")
w <- Variable(5, name = "w")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

slope <- Variable(1, name = "slope")
offset <- Variable(1, name = "offset")
quadratic_oceff <- Variable(1, name = "quadratic_coeff")

Tval <- 100
position <- Variable(2, Tval, name = "position")
velocity <- Variable(2, Tval, name = "velocity")
force <- Variable(2, Tval-1, name = "force")

xs <- Variable(80, name = "xs")
xsr <- Variable(200, name = "xsr")
xef <- Variable(80, name = "xef")

solvers <- installed_solvers()
solvers <- solvers[solvers %in% qp_solvers]
if("MOSEK" %in% installed_solvers())
  solvers <- c(solvers, "MOSEK")

solve_QP <- function(problem, solver_name) {
  solve(problem, solver = solver_name, verbose = TRUE)
}
